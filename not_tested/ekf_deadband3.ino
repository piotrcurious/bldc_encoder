 #include <Arduino.h>
#include <BasicLinearAlgebra.h>

// Pin definitions
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int COMMON_PIN = 9;

// Position tracking
volatile long estimated_mechanical_steps = 0;

// Sampling configuration
const int NUM_SAMPLES = 16;

// Control parameters
const float ANGLE_UPDATE_DEADBAND = 0.02f;     // ~1 degree
const float VELOCITY_STEP_THRESHOLD = 0.2f;    // rad/s
const float ZERO_VEL_LOCK_TIMEOUT = 500e3f;    // microseconds

// Excitation parameters
int EXCITATION_PWM_VALUE = 1;
const int EXCITATION_PULSE_WIDTH_US = 1000;
const int FLOATING_SETTLING_TIME_US = 10;

// Global filtered phase values
float g_valA = 0.0f;
float g_valB = 0.0f;
float g_valC = 0.0f;

// --- First Layer: Kalman Filter for ADC Smoothing ---
#define N_STATES_KF1 3

BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> F_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> H_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> Q_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> R_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> I_matrix_kf1;
BLA::Matrix<N_STATES_KF1, 1> z_kf1;

// KF1 tuning parameters
const float Q_DIAGONAL_VALUE_KF1 = 1.0f;
const float R_DIAGONAL_VALUE_KF1 = 100.0f;

// --- Second Layer: Extended Kalman Filter for Angle/Velocity ---
#define N_STATES_EKF 2

BLA::Matrix<N_STATES_EKF, 1> x_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_matrix_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> Q_ekf;
BLA::Matrix<1, 1> R_ekf;

// EKF tuning parameters
const float Q_ANGLE_EKF = 0.001f;
const float Q_VEL_EKF = 0.1f;
const float R_ANGLE_EKF = 0.5f;

// EKF state tracking
bool ekfLocked = false;
unsigned long lastNonZeroVelMicros = 0;
unsigned long prev_ekf_micros = 0;

// Auto-tuning parameters
const int TARGET_QUIESCENT_AVG_ADC = 80;
const int PWM_ADJUST_STEP = 1;
const int PWM_CALIBRATION_ITERATIONS = 2000;
const int PWM_CALIBRATION_TOLERANCE = 1;

/**
 * @brief Applies multivariate Kalman filter to phase measurements
 * @param raw_measurements Array of raw ADC readings [A, B, C]
 */
void applyMultivariateKalmanFilter(const int raw_measurements[]) {
    // Update measurement vector
    z_kf1(0) = static_cast<float>(raw_measurements[0]);
    z_kf1(1) = static_cast<float>(raw_measurements[1]);
    z_kf1(2) = static_cast<float>(raw_measurements[2]);

    // Prediction step
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_pred = P_kf1 + Q_kf1;
    
    // Update step
    BLA::Matrix<N_STATES_KF1, 1> y_k = z_kf1 - (H_kf1 * x_hat_kf1);
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S_k = H_kf1 * P_pred * ~H_kf1 + R_kf1;
    
    // Check for singularity
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S_k_inv = BLA::Inverse(S_k);
    if (isnan(S_k_inv(0,0)) || isinf(S_k_inv(0,0))) {
        return; // Skip update if singular
    }
    
    // Kalman gain and state update
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> K_k = P_pred * ~H_kf1 * S_k_inv;
    x_hat_kf1 = x_hat_kf1 + (K_k * y_k);
    P_kf1 = (I_matrix_kf1 - K_k * H_kf1) * P_pred;
}

/**
 * @brief Applies Extended Kalman Filter for angle and velocity estimation
 * @param inferred_angle_measurement Angle from atan2 calculation
 */
void applyEKF(float inferred_angle_measurement) {
    unsigned long current = micros();
    float dt = (current - prev_ekf_micros) / 1e6f;
    
    // Prevent division by zero and handle overflow
    if (dt <= 0 || dt > 1.0f) {
        prev_ekf_micros = current;
        return;
    }
    prev_ekf_micros = current;

    // Prediction step
    BLA::Matrix<2, 2> F = {1, dt, 0, 1};
    auto x_pred = F * x_ekf;
    auto P_pred = F * P_ekf * ~F + Q_ekf;

    // Calculate innovation with angle wrapping
    float ang_pred = x_pred(0);
    float y = inferred_angle_measurement - ang_pred;
    
    // Normalize angle difference to [-π, π]
    while (y > PI) y -= 2.0f * PI;
    while (y < -PI) y += 2.0f * PI;

    // Deadband and velocity-based update logic
    bool skip_update = (fabsf(y) < ANGLE_UPDATE_DEADBAND && 
                       fabsf(x_ekf(1)) < VELOCITY_STEP_THRESHOLD);
    
    if (skip_update) {
        // Handle zero-velocity locking
        if (!ekfLocked) {
            if (fabsf(x_ekf(1)) < VELOCITY_STEP_THRESHOLD) {
                if (current - lastNonZeroVelMicros > ZERO_VEL_LOCK_TIMEOUT) {
                    ekfLocked = true;
                }
            } else {
                lastNonZeroVelMicros = current;
            }
        }
        
        x_ekf = x_pred;
        P_ekf = P_pred;
        return;
    }

    // Reset lock state on significant motion
    ekfLocked = false;
    lastNonZeroVelMicros = current;

    // Update step
    BLA::Matrix<1, 2> H = {1, 0};
    auto S = H * P_pred * ~H + R_ekf;
    auto S_inv = BLA::Inverse(S);
    
    if (isnan(S_inv(0,0)) || isinf(S_inv(0,0))) {
        return; // Skip update if singular
    }
    
    auto K = P_pred * ~H * S_inv;
    BLA::Matrix<1, 1> yk = {y};
    x_ekf = x_pred + K * yk;
    P_ekf = (I_matrix_ekf - K * H) * P_pred;

    // Normalize angle to [0, 2π]
    while (x_ekf(0) >= 2.0f * PI) x_ekf(0) -= 2.0f * PI;
    while (x_ekf(0) < 0) x_ekf(0) += 2.0f * PI;
}

/**
 * @brief Performs reliable analog read with error checking
 * @param pin Analog pin to read
 * @return ADC reading (0-1023) or -1 on error
 */
int reliableAnalogRead(int pin) {
    if (pin < A0 || pin > A5) return -1; // Invalid pin
    return analogRead(pin);
}

/**
 * @brief Executes one sensing cycle: excite, float, measure
 */
void executeSensingCycle() {
    // Apply excitation
    analogWrite(COMMON_PIN, EXCITATION_PWM_VALUE);
    delayMicroseconds(EXCITATION_PULSE_WIDTH_US);

    // Switch to high impedance for BEMF sensing
    pinMode(COMMON_PIN, INPUT);
    delayMicroseconds(FLOATING_SETTLING_TIME_US);

    // Oversample for noise reduction
    long sumA = 0, sumB = 0, sumC = 0;
    int valid_samples = 0;
    
    for (int i = 0; i < NUM_SAMPLES; i++) {
        int readA = reliableAnalogRead(PHASE_A_PIN);
        int readB = reliableAnalogRead(PHASE_B_PIN);
        int readC = reliableAnalogRead(PHASE_C_PIN);
        
        if (readA >= 0 && readB >= 0 && readC >= 0) {
            sumA += readA;
            sumB += readB;
            sumC += readC;
            valid_samples++;
        }
    }

    // Handle case where no valid samples were obtained
    if (valid_samples == 0) {
        pinMode(COMMON_PIN, OUTPUT);
        digitalWrite(COMMON_PIN, LOW);
        return;
    }

    // Calculate averages
    int raw_vals[N_STATES_KF1];
    raw_vals[0] = sumA / valid_samples;
    raw_vals[1] = sumB / valid_samples;
    raw_vals[2] = sumC / valid_samples;

    // Apply first-layer Kalman filter
    applyMultivariateKalmanFilter(raw_vals);

    // Update global variables
    g_valA = x_hat_kf1(0);
    g_valB = x_hat_kf1(1);
    g_valC = x_hat_kf1(2);

    // Restore default state
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);
}

/**
 * @brief Calculates electrical angle from filtered phase voltages
 * @return Inferred electrical angle in radians
 */
float getInferredElectricalAngle() {
    // Clarke transformation components
    float x_component = (2.0f * g_valA - g_valB - g_valC);
    float y_component = (1.732050807568877f * (g_valB - g_valC)); // sqrt(3)

    // Handle degenerate case
    if (fabsf(x_component) < 1e-6f && fabsf(y_component) < 1e-6f) {
        return x_ekf(0); // Return current EKF estimate
    }

    return atan2(y_component, x_component);
}

/**
 * @brief Auto-tunes PWM excitation value for optimal sensing
 */
void autoTuneExcitationPwm() {
    Serial.println(F("\n--- Starting Auto-Tune ---"));
    Serial.println(F("Ensure motor is stationary!"));

    // Reset Kalman filter states
    x_hat_kf1.Fill(0.0f);
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        P_kf1(i, i) = 1000.0f; // High initial uncertainty
    }

    // Reset EKF states
    x_ekf.Fill(0.0f);
    P_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_EKF; i++) {
        P_ekf(i, i) = 10.0f; // High initial uncertainty
    }
    
    prev_ekf_micros = micros();
    ekfLocked = false;
    lastNonZeroVelMicros = micros();

    // Tuning loop
    for (int iteration = 0; iteration < PWM_CALIBRATION_ITERATIONS; iteration++) {
        executeSensingCycle();
        
        int current_avg_adc = static_cast<int>((g_valA + g_valB + g_valC) / 3.0f);
        int error = current_avg_adc - TARGET_QUIESCENT_AVG_ADC;

        if (iteration % 100 == 0) { // Reduce serial output frequency
            Serial.print(F("Step ")); Serial.print(iteration + 1);
            Serial.print(F(" PWM:")); Serial.print(EXCITATION_PWM_VALUE);
            Serial.print(F(" Avg:")); Serial.print(current_avg_adc);
            Serial.print(F(" Error:")); Serial.println(error);
        }

        if (abs(error) <= PWM_CALIBRATION_TOLERANCE) {
            Serial.println(F("--- Auto-tune complete! ---"));
            break;
        }

        // Adjust PWM value
        if (error < 0) {
            EXCITATION_PWM_VALUE += PWM_ADJUST_STEP;
        } else {
            EXCITATION_PWM_VALUE -= PWM_ADJUST_STEP;
        }

        // Clamp PWM value
        EXCITATION_PWM_VALUE = constrain(EXCITATION_PWM_VALUE, 0, 255);

        delay(10); // Reduced delay for faster tuning
    }

    Serial.print(F("Final PWM: ")); Serial.println(EXCITATION_PWM_VALUE);
    Serial.println(F("Resume operation - rotate motor for data"));
}

/**
 * @brief Initialize system matrices and parameters
 */
void initializeMatrices() {
    // Initialize identity matrices
    I_matrix_kf1.Fill(0.0f);
    I_matrix_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        I_matrix_kf1(i, i) = 1.0f;
    }
    for (int i = 0; i < N_STATES_EKF; i++) {
        I_matrix_ekf(i, i) = 1.0f;
    }

    // Initialize KF1 matrices
    F_kf1 = I_matrix_kf1; // Identity matrix
    H_kf1 = I_matrix_kf1; // Identity matrix
    
    Q_kf1.Fill(0.0f);
    R_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        Q_kf1(i, i) = Q_DIAGONAL_VALUE_KF1;
        R_kf1(i, i) = R_DIAGONAL_VALUE_KF1;
    }

    // Initialize EKF matrices
    Q_ekf.Fill(0.0f);
    Q_ekf(0, 0) = Q_ANGLE_EKF;
    Q_ekf(1, 1) = Q_VEL_EKF;
    
    R_ekf(0, 0) = R_ANGLE_EKF;
}

void setup() {
    Serial.begin(115200);
    while (!Serial) ; // Wait for serial port on Leonardo/Micro
    
    Serial.println(F("BLDC Encoder Initialized"));

    // Initialize matrices
    initializeMatrices();

    // Configure ADC for faster, more accurate readings
    analogReference(INTERNAL); // Use internal 1.1V reference
    ADCSRA &= ~((1 << ADPS0) | (1 << ADPS1) | (1 << ADPS2)); // Clear prescaler
    ADCSRA |= (1 << ADPS2) | (1 << ADPS0); // Set prescaler to 32 (500kHz)
    ADCSRA |= (1 << ADEN); // Enable ADC

    // Increase PWM frequency for smoother excitation
    TCCR1B = (TCCR1B & 0xF8) | 0x01; // Set prescaler to 1 for higher frequency

    // Initialize common pin
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);

    // Perform auto-tuning
    autoTuneExcitationPwm();

    // Initialize states after tuning
    x_hat_kf1.Fill(static_cast<float>(TARGET_QUIESCENT_AVG_ADC));
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        P_kf1(i, i) = 100.0f; // Lower uncertainty after tuning
    }

    // Initialize EKF with first measurement
    executeSensingCycle();
    x_ekf(0) = getInferredElectricalAngle();
    x_ekf(1) = 0.0f;
    P_ekf(0, 0) = 1.0f;
    P_ekf(1, 1) = 1.0f;
    
    prev_ekf_micros = micros();
    
    Serial.println(F("System ready"));
}

void loop() {
    // Main sensing and estimation cycle
    executeSensingCycle();
    float inferred_angle = getInferredElectricalAngle();
    applyEKF(inferred_angle);

    // Update mechanical position estimate based on velocity
    float velocity = x_ekf(1);
    if (fabsf(velocity) > VELOCITY_STEP_THRESHOLD) {
        estimated_mechanical_steps = round(x_ekf(0) / (2.0f * PI) * 6.0f);
    }

    // Determine rotation direction
    int rotation_direction = (velocity > VELOCITY_STEP_THRESHOLD) ? 1 :
                           (velocity < -VELOCITY_STEP_THRESHOLD) ? -1 : 0;

    // Output data for Serial Plotter (tab-separated)
    Serial.print(g_valA, 2); Serial.print(F("\t"));
    Serial.print(g_valB, 2); Serial.print(F("\t"));
    Serial.print(g_valC, 2); Serial.print(F("\t"));
    Serial.print(inferred_angle * 180.0f / PI, 2); Serial.print(F("\t"));
    Serial.print(x_ekf(0) * 180.0f / PI, 2); Serial.print(F("\t"));
    Serial.print(x_ekf(1), 2); Serial.print(F("\t"));
    Serial.print(estimated_mechanical_steps); Serial.print(F("\t"));
    Serial.println(rotation_direction);

    // Handle serial commands
    if (Serial.available()) {
        char command = Serial.read();
        Serial.read(); // Clear any remaining characters
        
        if (command == 'T' || command == 't') {
            Serial.println(F("Re-tuning..."));
            autoTuneExcitationPwm();
            
            // Reset states after re-tuning
            x_ekf(0) = getInferredElectricalAngle();
            x_ekf(1) = 0.0f;
            P_ekf(0, 0) = 1.0f;
            P_ekf(1, 1) = 1.0f;
            prev_ekf_micros = micros();
            estimated_mechanical_steps = 0;
            ekfLocked = false;
        }
    }
}
