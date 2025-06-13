#include <Arduino.h>
#include <BasicLinearAlgebra.h>

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~ Configuration Namespace ~~
// ~~ All user-tunable parameters are located here for easy access. ~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace Config {
    // --- Pin Definitions ---
    constexpr int PHASE_A_PIN = A0;
    constexpr int PHASE_B_PIN = A1;
    constexpr int PHASE_C_PIN = A2;
    constexpr int COMMON_PIN = 9; // Must be a PWM pin (Timer1 for this config: 9 or 10 on Uno)

    // --- Physical & Mathematical Constants ---
    constexpr float ELECTRICAL_STEPS_PER_REV = 6.0f; // Steps per electrical revolution (motor specific)
   constexpr float TWO_PI_F = 2.0f * PI;
    constexpr float SQRT_3 = 1.73205080757f;

    // --- Sampling Configuration ---
    constexpr int NUM_SAMPLES = 16; // Number of ADC samples to average

    // --- EKF Control Parameters ---
    constexpr float ANGLE_UPDATE_DEADBAND = 0.02f;      // rad (~1 degree), for noise rejection at standstill
    constexpr float VELOCITY_STEP_THRESHOLD = 0.2f;     // rad/s, minimum velocity to be considered "moving"
    constexpr unsigned long ZERO_VEL_LOCK_TIMEOUT = 500000UL; // microseconds

    // --- Excitation Parameters ---
    int EXCITATION_PWM_VALUE = 1; // Initial value, will be auto-tuned
    constexpr int EXCITATION_PULSE_WIDTH_US = 1000;
    constexpr int FLOATING_SETTLING_TIME_US = 10;

    // --- Auto-Tuning Parameters ---
    constexpr int TARGET_QUIESCENT_AVG_ADC = 80;    // Target ADC value for phase sense when motor is still
    constexpr int PWM_ADJUST_STEP = 1;              // PWM adjustment step during tuning
    constexpr int PWM_CALIBRATION_ITERATIONS = 2000;
    constexpr int PWM_CALIBRATION_TOLERANCE = 1;

    // --- KF1 (ADC Smoothing) Tuning ---
    constexpr float Q_DIAGONAL_VALUE_KF1 = 1.0f;    // Process noise - higher means more trust in measurements
    constexpr float R_DIAGONAL_VALUE_KF1 = 100.0f;  // Measurement noise - higher means more trust in model

    // --- EKF (Angle/Velocity) Tuning ---
    constexpr float Q_ANGLE_EKF = 0.001f;           // Process noise for angle
    constexpr float Q_VEL_EKF = 0.1f;               // Process noise for velocity
    constexpr float R_ANGLE_EKF = 0.5f;             // Measurement noise for angle
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~ Global State & Matrices ~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// --- Layer 1: Kalman Filter for ADC Smoothing ---
#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1;

// --- Layer 2: Extended Kalman Filter for Angle/Velocity ---
#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf;
BLA::Matrix<1, 1> R_ekf;

// --- EKF State & Position Tracking ---
bool ekfLocked = false;
unsigned long lastNonZeroVelMicros = 0;
unsigned long prev_ekf_micros = 0;
long revolution_count = 0;
float unwrapped_electrical_angle = 0.0f;
long estimated_mechanical_steps = 0;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~ Function Prototypes ~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void initializeMatrices();
void autoTuneExcitationPwm();
void executeSensingCycle();
void applyMultivariateKalmanFilter(const int raw_measurements[]);
float getInferredElectricalAngle();
void applyEKF(float inferred_angle_measurement);
void logData();
int checkedAnalogRead(int pin);
float normalizeAngle(float angle);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~ Main Setup & Loop ~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void setup() {
    Serial.begin(115200);
    while (!Serial) { delay(10); } // Wait for serial port on some Arduinos
    
    Serial.println(F("Initializing BLDC Sensorless Encoder..."));

    initializeMatrices();

    // --- Hardware Configuration ---
    // Use the stable internal 1.1V reference. Note: Max input voltage on ADC pins is now 1.1V!
    analogReference(INTERNAL);
    
    // Set ADC prescaler to 32 for a 500kHz ADC clock (16MHz / 32).
    // Note: This is faster than the datasheet's recommended 200kHz max, trading some
    // precision for speed, which is acceptable for this differential application.
    ADCSRA &= ~((1 << ADPS2) | (1 << ADPS1) | (1 << ADPS0));
    ADCSRA |= (1 << ADPS2) | (1 << ADPS0);

    // Set Timer1 prescaler to 1 for a higher PWM frequency on pins 9 & 10 (~31.25 kHz).
    // This makes the excitation pulse quieter and potentially more effective.
    TCCR1B = (TCCR1B & 0b11111000) | 0x01;

    pinMode(Config::COMMON_PIN, OUTPUT);
    digitalWrite(Config::COMMON_PIN, LOW);

    autoTuneExcitationPwm();

    // Initialize states after tuning
    x_hat_kf1.Fill(static_cast<float>(Config::TARGET_QUIESCENT_AVG_ADC));
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) { P_kf1(i, i) = 100.0f; }

    // Initialize EKF with the first valid measurement
    executeSensingCycle();
    x_ekf(0) = getInferredElectricalAngle();
    x_ekf(1) = 0.0f;
    P_ekf.Fill(0.0f);
    P_ekf(0, 0) = 1.0f;
    P_ekf(1, 1) = 1.0f;
    
    prev_ekf_micros = micros();
    unwrapped_electrical_angle = x_ekf(0); // Initialize unwrapped angle
    
    Serial.println(F("System ready. Rotate motor to begin estimation."));
}

void loop() {
    // Main sensing and estimation cycle
    executeSensingCycle();
    float inferred_angle = getInferredElectricalAngle();
    applyEKF(inferred_angle);

    // Update mechanical position estimate from the continuous unwrapped angle
    estimated_mechanical_steps = round(unwrapped_electrical_angle / Config::TWO_PI_F * Config::ELECTRICAL_STEPS_PER_REV);
    
    logData();

    // Handle serial commands for re-tuning
    if (Serial.available() > 0) {
        char command = toupper(Serial.read());
        // Clear remaining characters in the serial buffer
        while (Serial.available() > 0) { Serial.read(); }
        
        if (command == 'T') {
            Serial.println(F("\nRe-tuning requested..."));
            autoTuneExcitationPwm();
            
            // Reset states after re-tuning
            executeSensingCycle();
            x_ekf(0) = getInferredElectricalAngle();
            x_ekf(1) = 0.0f;
            P_ekf(0, 0) = 1.0f;
            P_ekf(1, 1) = 1.0f;
            prev_ekf_micros = micros();
            revolution_count = 0;
            unwrapped_electrical_angle = x_ekf(0);
            ekfLocked = false;
            Serial.println(F("Re-tuning complete. System ready."));
        }
    }
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~ Function Implementations ~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
 * @brief Initializes KF and EKF matrices with configured values.
 */
void initializeMatrices() {
    I_matrix_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) I_matrix_kf1(i, i) = 1.0f;
    
    I_matrix_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_EKF; i++) I_matrix_ekf(i, i) = 1.0f;

    Q_kf1.Fill(0.0f);
    R_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        Q_kf1(i, i) = Config::Q_DIAGONAL_VALUE_KF1;
        R_kf1(i, i) = Config::R_DIAGONAL_VALUE_KF1;
    }

    Q_ekf.Fill(0.0f);
    Q_ekf(0, 0) = Config::Q_ANGLE_EKF;
    Q_ekf(1, 1) = Config::Q_VEL_EKF;
    
    R_ekf(0, 0) = Config::R_ANGLE_EKF;
}


/**
 * @brief Auto-tunes the PWM excitation value to achieve the target quiescent ADC reading.
 */
void autoTuneExcitationPwm() {
    Serial.println(F("\n--- Starting Auto-Tune ---"));
    Serial.println(F("Ensure motor is stationary!"));

    // Reset filters to a high-uncertainty state
    x_hat_kf1.Fill(0.0f);
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) P_kf1(i, i) = 1000.0f;

    for (int iteration = 0; iteration < Config::PWM_CALIBRATION_ITERATIONS; iteration++) {
        executeSensingCycle();
        
        int current_avg_adc = static_cast<int>((x_hat_kf1(0) + x_hat_kf1(1) + x_hat_kf1(2)) / 3.0f);
        int error = current_avg_adc - Config::TARGET_QUIESCENT_AVG_ADC;

        if (iteration % 100 == 0) { // Reduce serial output frequency
            Serial.print(F("Step ")); Serial.print(iteration);
            Serial.print(F(" | PWM: ")); Serial.print(Config::EXCITATION_PWM_VALUE);
            Serial.print(F(" | Avg ADC: ")); Serial.print(current_avg_adc);
            Serial.print(F(" | Error: ")); Serial.println(error);
        }

        if (abs(error) <= Config::PWM_CALIBRATION_TOLERANCE) {
            Serial.println(F("--- Auto-tune complete! ---"));
            break;
        }

        Config::EXCITATION_PWM_VALUE += (error < 0) ? Config::PWM_ADJUST_STEP : -Config::PWM_ADJUST_STEP;
        Config::EXCITATION_PWM_VALUE = constrain(Config::EXCITATION_PWM_VALUE, 1, 255); // PWM must be > 0

        delay(5); // Short delay for stability
    }
    Serial.print(F("Final Tuned PWM Value: ")); Serial.println(Config::EXCITATION_PWM_VALUE);
}

/**
 * @brief Executes one sensing cycle: excites phases, waits, and measures response.
 * Updates the first-layer Kalman filter state (x_hat_kf1).
 */
void executeSensingCycle() {
    analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
    delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);

    pinMode(Config::COMMON_PIN, INPUT); // Set pin to high-Z to sense BEMF
    delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);

    long sumA = 0, sumB = 0, sumC = 0;
    int valid_samples = 0;
    for (int i = 0; i < Config::NUM_SAMPLES; i++) {
        int readA = checkedAnalogRead(Config::PHASE_A_PIN);
        int readB = checkedAnalogRead(Config::PHASE_B_PIN);
        int readC = checkedAnalogRead(Config::PHASE_C_PIN);
        
        if (readA >= 0 && readB >= 0 && readC >= 0) {
            sumA += readA; sumB += readB; sumC += readC;
            valid_samples++;
        }
    }
    
    pinMode(Config::COMMON_PIN, OUTPUT); // Restore pin to output
    digitalWrite(Config::COMMON_PIN, LOW);

    if (valid_samples > 0) {
        int raw_vals[N_STATES_KF1] = {
            static_cast<int>(sumA / valid_samples),
            static_cast<int>(sumB / valid_samples),
            static_cast<int>(sumC / valid_samples)
        };
        applyMultivariateKalmanFilter(raw_vals);
    }
}

/**
 * @brief Applies a multivariate Kalman filter to smooth raw ADC phase measurements.
 * @param raw_measurements Array of averaged ADC readings [A, B, C].
 */
void applyMultivariateKalmanFilter(const int raw_measurements[]) {
    BLA::Matrix<N_STATES_KF1, 1> z_kf1 = {
        static_cast<float>(raw_measurements[0]),
        static_cast<float>(raw_measurements[1]),
        static_cast<float>(raw_measurements[2])
    };

    // Prediction (state is assumed static between measurements, F=I)
    auto P_pred = P_kf1 + Q_kf1;

    // Update
    auto y_k = z_kf1 - x_hat_kf1; // H=I
    auto S_k = P_pred + R_kf1;    // H=I
    auto S_k_inv = BLA::Inverse(S_k);

    if (isnan(S_k_inv(0, 0)) || isinf(S_k_inv(0, 0))) { return; } // Skip update if singular

    auto K_k = P_pred * S_k_inv;
    x_hat_kf1 = x_hat_kf1 + K_k * y_k;
    P_kf1 = (I_matrix_kf1 - K_k) * P_pred;
}

/**
 * @brief Infers the electrical angle from filtered phase voltages using Clarke Transform.
 * @return Inferred electrical angle in radians [0, 2PI].
 */
float getInferredElectricalAngle() {
    float valA = x_hat_kf1(0), valB = x_hat_kf1(1), valC = x_hat_kf1(2);

    float alpha = (2.0f * valA - valB - valC) / 3.0f;
    float beta = (Config::SQRT_3 * (valB - valC)) / 3.0f;

    // Return current EKF angle if magnitude is too small (prevents instability)
    if (fabsf(alpha) < 1e-6f && fabsf(beta) < 1e-6f) {
        return x_ekf(0);
    }

    float angle = atan2(beta, alpha);
    return (angle < 0.0f) ? angle + Config::TWO_PI_F : angle; // Map to [0, 2PI]
}

/**
 * @brief Applies an Extended Kalman Filter to estimate angle and velocity.
 * @param inferred_angle_measurement The raw angle from getInferredElectricalAngle().
 */
void applyEKF(float inferred_angle_measurement) {
    unsigned long current_micros = micros();
    float dt = (current_micros - prev_ekf_micros) / 1e6f;
    
    if (dt <= 0.0f || dt > 0.5f) { // Ignore spurious timings
        prev_ekf_micros = current_micros;
        return;
    }
    prev_ekf_micros = current_micros;

    float last_angle = x_ekf(0); // Store angle from previous state for unwrapping

    // Prediction
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> F = {1, dt, 0, 1};
    auto x_pred = F * x_ekf;
    auto P_pred = F * P_ekf * ~F + Q_ekf;

    // Calculate innovation (error) with angle wrapping
    float y = inferred_angle_measurement - x_pred(0);
    y = normalizeAngle(y); // Normalize innovation to [-PI, PI]

    // Skip update if measurement is too close to prediction (prevents drift at standstill)
    bool skip_update = (fabsf(y) < Config::ANGLE_UPDATE_DEADBAND && fabsf(x_ekf(1)) < Config::VELOCITY_STEP_THRESHOLD);

    if (skip_update) {
        if (!ekfLocked) {
            if (fabsf(x_ekf(1)) < Config::VELOCITY_STEP_THRESHOLD) {
                if (current_micros - lastNonZeroVelMicros > Config::ZERO_VEL_LOCK_TIMEOUT) {
                    ekfLocked = true;
                }
            } else {
                lastNonZeroVelMicros = current_micros;
            }
        }
        x_ekf = x_pred;
        P_ekf = P_pred;
        return;
    }

    // Reset lock state on significant motion
    ekfLocked = false;
    lastNonZeroVelMicros = current_micros;

    // Update
    BLA::Matrix<1, N_STATES_EKF> H = {1, 0};
    auto S = H * P_pred * ~H + R_ekf;
    auto S_inv = BLA::Inverse(S);
    if (isnan(S_inv(0, 0)) || isinf(S_inv(0, 0))) { return; } // Skip if singular

    auto K = P_pred * ~H * S_inv;
    BLA::Matrix<1,1> y_mat = {y};

    x_ekf = x_pred + K * y_mat;
    P_ekf = (I_matrix_ekf - K * H) * P_pred;

    // Normalize final angle estimate to [0, 2PI]
    x_ekf(0) = fmodf(x_ekf(0), Config::TWO_PI_F);
    if (x_ekf(0) < 0.0f) x_ekf(0) += Config::TWO_PI_F;

    // --- Angle Unwrapping for Continuous Position Tracking ---
    float angle_diff = x_ekf(0) - last_angle;
    if (angle_diff > PI) {
        revolution_count--; // Wrapped from near 2PI to 0
    } else if (angle_diff < -PI) {
        revolution_count++; // Wrapped from near 0 to 2PI
    }
    unwrapped_electrical_angle = x_ekf(0) + revolution_count * Config::TWO_PI_F;
}

/**
 * @brief Prints system state data to the Serial Plotter.
 */
void logData() {
    static unsigned long last_log_time = 0;
    if (millis() - last_log_time < 20) return; // Limit output rate
    last_log_time = millis();

    float velocity = x_ekf(1);
    int direction = (velocity > Config::VELOCITY_STEP_THRESHOLD) ? 1 : (velocity < -Config::VELOCITY_STEP_THRESHOLD) ? -1 : 0;

    // Format: PhaseA, PhaseB, PhaseC, InferredAngle, EKFAngle, EKFVelocity, MechSteps, Direction
    Serial.print(x_hat_kf1(0), 2); Serial.print(F("\t"));
    Serial.print(x_hat_kf1(1), 2); Serial.print(F("\t"));
    Serial.print(x_hat_kf1(2), 2); Serial.print(F("\t"));
    Serial.print(getInferredElectricalAngle() * 180.0f / PI, 2); Serial.print(F("\t"));
    Serial.print(x_ekf(0) * 180.0f / PI, 2); Serial.print(F("\t"));
    Serial.print(velocity, 2); Serial.print(F("\t"));
    Serial.print(estimated_mechanical_steps); Serial.print(F("\t"));
    Serial.println(direction);
}


/**
 * @brief Reads an analog pin after checking its validity for the board.
 * @param pin The analog pin to read (e.g., A0, A1).
 * @return The ADC value (0-1023) or -1 if the pin is invalid.
 */
int checkedAnalogRead(int pin) {
    // Note: This check is for ATmega328p-based boards (Uno, Nano).
    // It may need adjustment for other boards like the Mega.
    if (pin >= 14 && pin <= 19) {
        return analogRead(pin);
    }
    return -1;
}

/**
 * @brief Normalizes an angle to the range [-PI, PI].
 * @param angle The angle in radians.
 * @return The normalized angle in radians.
 */
float normalizeAngle(float angle) {
    angle = fmodf(angle + PI, Config::TWO_PI_F);
    if (angle < 0.0f) {
        angle += Config::TWO_PI_F;
    }
    return angle - PI;
}
