#include <Arduino.h>
#include <BasicLinearAlgebra.h> // Include the BLA library

const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int COMMON_PIN = 9; // Pin 9 is controlled by Timer1 (along with pin 10)

// Renamed from encoderPosition, now represents the estimated mechanical angle in 'steps'
// Each full electrical rotation (2*PI radians) corresponds to 6 'Hall steps'
volatile long estimated_mechanical_steps = 0; // Continuous estimate of mechanical angle in Hall steps

const int NUM_SAMPLES = 16;

// Excitation parameters
int EXCITATION_PWM_VALUE = 1; // THIS WILL BE AUTO-TUNED. Initial guess for PWM strength (0-255).
const int EXCITATION_PULSE_WIDTH_US = 1000; // Duration of the initial PWM pulse on COMMON_PIN.
const int FLOATING_SETTLING_TIME_US = 10;     // Crucial: Time to float after excitation pulse before sensing.

// Global variables for plotting (current averaged ADC values from first KF)
float g_valA = 0.0f;
float g_valB = 0.0f;
float g_valC = 0.0f;

// --- Velocity Smoothing and Spurious Rotation Prevention ---
const int VELOCITY_HISTORY_SIZE = 10;
float velocity_history[VELOCITY_HISTORY_SIZE];
int velocity_history_index = 0;
bool velocity_history_filled = false;

// Improved thresholds for rotation detection
const float ROTATION_VELOCITY_THRESHOLD = 0.5f;  // Increased from 0.1f rad/s
const float STATIONARY_VELOCITY_THRESHOLD = 0.2f; // Hysteresis threshold for stopping
const int MIN_CONSISTENT_SAMPLES = 5; // Minimum samples showing consistent direction
const float MAX_STATIONARY_VELOCITY_STD = 0.1f; // Max std dev when stationary

// Rotation state tracking
enum RotationState {
    STATIONARY,
    ROTATING_CW,
    ROTATING_CCW
};

RotationState current_rotation_state = STATIONARY;
int consistent_direction_count = 0;
float last_stable_angle = 0.0f;

// --- First Layer: Kalman Filter Variables for ADC Smoothing ---
#define N_STATES_KF1 3 // Number of states (V_A, V_B, V_C)

BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1; // State vector [V_A, V_B, V_C]'
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> F_kf1; // State Transition Matrix - Identity
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> H_kf1; // Observation Matrix - Identity

const float Q_DIAGONAL_VALUE_KF1 = 1.0f; // Process Noise Covariance for first KF
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> Q_kf1;

const float R_DIAGONAL_VALUE_KF1 = 100.0f; // Measurement Noise Covariance for first KF
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> R_kf1;

BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1; // Initial Estimate Covariance Matrix for first KF

// Explicit Identity matrix for BLA operations - now properly used
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> I_matrix_kf1;

// Measurement vector (z_k) for first KF
BLA::Matrix<N_STATES_KF1, 1> z_kf1;

/**
 * @brief Applies a multivariate Kalman filter to the phase measurements.
 * Updates the global x_hat_kf1 (state estimate) and P_kf1 (covariance) matrices.
 * @param raw_measurements An array of raw ADC readings [A, B, C].
 */
void applyMultivariateKalmanFilter(int raw_measurements[]) {
    z_kf1(0) = (float)raw_measurements[0];
    z_kf1(1) = (float)raw_measurements[1];
    z_kf1(2) = (float)raw_measurements[2];

    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_pred = P_kf1 + Q_kf1;

    BLA::Matrix<N_STATES_KF1, 1> y_k = z_kf1 - (H_kf1 * x_hat_kf1);
    // Corrected Transpose syntax
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S_k = H_kf1 * P_pred * ~H_kf1 + R_kf1; // Use ~ operator for transpose

    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S_k_inv = BLA::Inverse(S_k);
    if (isnan(S_k_inv(0,0))) { return; } // Skip update if S_k is singular

    // Corrected Transpose syntax
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> K_k = P_pred * ~H_kf1 * S_k_inv; // Use ~ operator for transpose

    x_hat_kf1 = x_hat_kf1 + (K_k * y_k);
    // Corrected Identity syntax: Use the pre-initialized I_matrix_kf1
    P_kf1 = (I_matrix_kf1 - K_k * H_kf1) * P_pred; // Updated covariance
}

// --- Second Layer: Extended Kalman Filter (EKF) for Angle and Velocity ---
#define N_STATES_EKF 2 // States: [electrical_angle, electrical_velocity]

BLA::Matrix<N_STATES_EKF, 1> x_ekf; // EKF State vector [electrical_angle, electrical_velocity]'
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf; // EKF Covariance matrix

// Explicit Identity matrix for BLA operations for EKF
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_matrix_ekf;

// Process Noise Covariance Matrix (Q_ekf) - REDUCED for better stability
// Smaller values make the filter less responsive but more stable
const float Q_ANGLE_EKF = 0.0001f; // Reduced from 0.001f for more stable angle
const float Q_VEL_EKF = 0.01f;     // Reduced from 0.1f for less noisy velocity
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> Q_ekf = {
    Q_ANGLE_EKF, 0,
    0, Q_VEL_EKF
};

// Measurement Noise Covariance Matrix (R_ekf) for the inferred_angle measurement
const float R_ANGLE_EKF = 0.05f; // Reduced from 0.1f for trusting measurements more
BLA::Matrix<1, 1> R_ekf = {R_ANGLE_EKF}; // Measurement is a single value

unsigned long prev_ekf_micros = 0; // For calculating dt

/**
 * @brief Adds a velocity sample to the history buffer and checks for spurious rotation.
 * @param velocity The current velocity estimate from EKF.
 */
void updateVelocityHistory(float velocity) {
    velocity_history[velocity_history_index] = velocity;
    velocity_history_index = (velocity_history_index + 1) % VELOCITY_HISTORY_SIZE;
    
    if (!velocity_history_filled && velocity_history_index == 0) {
        velocity_history_filled = true;
    }
}

/**
 * @brief Calculates the average velocity from the history buffer.
 * @return Average velocity over the history window.
 */
float getAverageVelocity() {
    if (!velocity_history_filled && velocity_history_index == 0) {
        return 0.0f; // No data yet
    }
    
    float sum = 0.0f;
    int count = velocity_history_filled ? VELOCITY_HISTORY_SIZE : velocity_history_index;
    
    for (int i = 0; i < count; i++) {
        sum += velocity_history[i];
    }
    
    return sum / count;
}

/**
 * @brief Calculates the standard deviation of velocity from the history buffer.
 * @return Standard deviation of velocity over the history window.
 */
float getVelocityStandardDeviation() {
    if (!velocity_history_filled && velocity_history_index < 2) {
        return 0.0f; // Need at least 2 samples
    }
    
    float avg = getAverageVelocity();
    float sum_sq_diff = 0.0f;
    int count = velocity_history_filled ? VELOCITY_HISTORY_SIZE : velocity_history_index;
    
    for (int i = 0; i < count; i++) {
        float diff = velocity_history[i] - avg;
        sum_sq_diff += diff * diff;
    }
    
    return sqrt(sum_sq_diff / count);
}

/**
 * @brief Applies the Extended Kalman Filter to estimate electrical angle and velocity.
 * @param inferred_angle_measurement The angle derived from g_valA, g_valB, g_valC using atan2.
 */
void applyEKF(float inferred_angle_measurement) {
    unsigned long current_micros = micros();
    float dt = (current_micros - prev_ekf_micros) / 1000000.0f; // Convert to seconds
    prev_ekf_micros = current_micros;

    if (dt == 0 || dt > 0.1f) return; // Avoid division by zero or excessive dt

    // --- EKF Prediction Step ---
    // State Transition Matrix (F_ekf) - Jacobian of f(x)
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> F_ekf = {
        1, dt,  // angle_k = angle_k-1 + velocity_k-1 * dt
        0, 1    // velocity_k = velocity_k-1
    };

    BLA::Matrix<N_STATES_EKF, 1> x_pred = F_ekf * x_ekf;
    // Corrected Transpose syntax
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_pred = F_ekf * P_ekf * ~F_ekf + Q_ekf; // Use ~ operator for transpose

    // --- EKF Update Step ---
    // Measurement (z_ekf)
    BLA::Matrix<1, 1> z_ekf = {inferred_angle_measurement};

    // Predicted Measurement (h(x_pred))
    // Here, h(x) = x[0] (the angle itself)
    BLA::Matrix<1, 1> h_x = {x_pred(0)};

    // Jacobian of the Measurement Model (H_ekf)
    BLA::Matrix<1, N_STATES_EKF> H_ekf = {1, 0}; // dh/d(angle) = 1, dh/d(velocity) = 0

    BLA::Matrix<1, 1> y = z_ekf - h_x; // Innovation (measurement residual)

    // Handle angle wrapping for innovation: Ensure y is between -PI and PI
    while (y(0) > PI) y(0) -= 2 * PI;
    while (y(0) < -PI) y(0) += 2 * PI;

    // Corrected Transpose syntax
    BLA::Matrix<1, 1> S = H_ekf * P_pred * ~H_ekf + R_ekf; // Innovation covariance // Use ~ operator for transpose
    BLA::Matrix<1, 1> S_inv = BLA::Inverse(S);
    if (isnan(S_inv(0,0))) { return; } // Skip update if S is singular

    // Corrected Transpose syntax
    BLA::Matrix<N_STATES_EKF, 1> K = P_pred * ~H_ekf * S_inv; // Kalman Gain // Use ~ operator for transpose

    x_ekf = x_pred + K * y; // Updated state
    // Corrected Identity syntax: Use the pre-initialized I_matrix_ekf
    P_ekf = (I_matrix_ekf - K * H_ekf) * P_pred; // Updated covariance

    // Ensure angle wraps correctly in the state vector (0 to 2*PI)
    while (x_ekf(0) >= 2 * PI) x_ekf(0) -= 2 * PI;
    while (x_ekf(0) < 0) x_ekf(0) += 2 * PI;

    // Update velocity history for spurious rotation detection
    updateVelocityHistory(x_ekf(1));
}

/**
 * @brief Determines rotation state with improved spurious rotation prevention.
 * @return Rotation direction: -1 (CCW), 0 (Stationary), 1 (CW)
 */
int getRotationDirection() {
    float current_velocity = x_ekf(1);
    float avg_velocity = getAverageVelocity();
    float velocity_std = getVelocityStandardDeviation();
    
    // Check if we're in a stationary state based on standard deviation
    bool likely_stationary = (velocity_std < MAX_STATIONARY_VELOCITY_STD);
    
    // State machine for rotation detection with hysteresis
    switch (current_rotation_state) {
        case STATIONARY:
            // Only transition out of stationary if we have consistent direction
            if (fabs(avg_velocity) > ROTATION_VELOCITY_THRESHOLD && !likely_stationary) {
                if (avg_velocity > 0) {
                    consistent_direction_count++;
                    if (consistent_direction_count >= MIN_CONSISTENT_SAMPLES) {
                        current_rotation_state = ROTATING_CW;
                        consistent_direction_count = 0;
                        return 1;
                    }
                } else {
                    consistent_direction_count++;
                    if (consistent_direction_count >= MIN_CONSISTENT_SAMPLES) {
                        current_rotation_state = ROTATING_CCW;
                        consistent_direction_count = 0;
                        return -1;
                    }
                }
            } else {
                consistent_direction_count = 0; // Reset if not consistent
            }
            return 0;
            
        case ROTATING_CW:
            // Stay in CW state unless velocity drops below lower threshold or changes direction
            if (fabs(avg_velocity) < STATIONARY_VELOCITY_THRESHOLD || likely_stationary) {
                current_rotation_state = STATIONARY;
                consistent_direction_count = 0;
                return 0;
            } else if (avg_velocity < -ROTATION_VELOCITY_THRESHOLD) {
                consistent_direction_count++;
                if (consistent_direction_count >= MIN_CONSISTENT_SAMPLES) {
                    current_rotation_state = ROTATING_CCW;
                    consistent_direction_count = 0;
                    return -1;
                }
            } else {
                consistent_direction_count = 0;
            }
            return 1;
            
        case ROTATING_CCW:
            // Stay in CCW state unless velocity drops below lower threshold or changes direction
            if (fabs(avg_velocity) < STATIONARY_VELOCITY_THRESHOLD || likely_stationary) {
                current_rotation_state = STATIONARY;
                consistent_direction_count = 0;
                return 0;
            } else if (avg_velocity > ROTATION_VELOCITY_THRESHOLD) {
                consistent_direction_count++;
                if (consistent_direction_count >= MIN_CONSISTENT_SAMPLES) {
                    current_rotation_state = ROTATING_CW;
                    consistent_direction_count = 0;
                    return 1;
                }
            } else {
                consistent_direction_count = 0;
            }
            return -1;
    }
    
    return 0; // Default to stationary
}

// Auto-tuning parameters for the quiescent phase average
const int TARGET_QUIESCENT_AVG_ADC = 80; // Target for average phase voltage in stationary state (0-1023 scale)
const int PWM_ADJUST_STEP = 1 ;            // Fine-grained adjustment for PWM value
const int PWM_CALIBRATION_ITERATIONS = 2000; // Increased iterations for better convergence
const int PWM_CALIBRATION_TOLERANCE = 1;  // Tighter acceptable deviation for tuning

/**
 * @brief Performs a single analog read.
 * @param pin The analog pin to read (A0, A1, A2).
 * @return The 10-bit ADC reading (0-1023).
 */
int reliableAnalogRead(int pin) {
    return analogRead(pin);
}

/**
 * @brief Performs a single cycle of PWM excitation, floating, and sensing.
 * This function calculates the average phase voltages (g_valA/B/C) but does not
 * determine the Hall state.
 */
void executeSensingCycle() {
    // 1. Apply PWM excitation to COMMON_PIN
    analogWrite(COMMON_PIN, EXCITATION_PWM_VALUE);
    delayMicroseconds(EXCITATION_PULSE_WIDTH_US);

    // 2. Immediately set COMMON_PIN to INPUT for BEMF sensing
    pinMode(COMMON_PIN, INPUT); // Set to high impedance for sensing

    delayMicroseconds(FLOATING_SETTLING_TIME_US);

    // 3. Aggressive oversampling with reliable analog reads
    long sumA = 0, sumB = 0, sumC = 0;
    for (int i = 0; i < NUM_SAMPLES; i++) {
        sumA += reliableAnalogRead(PHASE_A_PIN);
        sumB += reliableAnalogRead(PHASE_B_PIN);
        sumC += reliableAnalogRead(PHASE_C_PIN);
    }

    // Get raw averaged values
    int raw_vals[N_STATES_KF1];
    raw_vals[0] = sumA / NUM_SAMPLES;
    raw_vals[1] = sumB / NUM_SAMPLES;
    raw_vals[2] = sumC / NUM_SAMPLES;

    // Apply First Layer Kalman Filter
    applyMultivariateKalmanFilter(raw_vals);

    // Update global variables with filtered estimates
    g_valA = x_hat_kf1(0);
    g_valB = x_hat_kf1(1);
    g_valC = x_hat_kf1(2);

    // 4. Re-establish default state: COMMON_PIN as OUTPUT, LOW (stable ground sink)
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);
}

/**
 * @brief Calculates an inferred electrical angle from the filtered phase voltages.
 * This acts as the measurement for the second-layer EKF.
 * @return Inferred electrical angle in radians (-PI to PI).
 */
float getInferredElectricalAngle() {
    // This formula assumes a sinusoidal relationship between phases and angle.
    // It's a common approximation for sensorless BLDC position estimation.
    // X-component: (2*Va - Vb - Vc)
    // Y-component: sqrt(3) * (Vb - Vc)
    // Adjust phases based on your motor's specific winding and sensor setup if needed.
    // For typical BLDC, phases are 120 deg apart. This maps them to a rotating vector.
    float x_component = (2.0f * g_valA - g_valB - g_valC);
    float y_component = (1.732f * (g_valB - g_valC)); // SQRT3 is ~1.732

    // Avoid atan2(0,0) if all values are identical or zero, though unlikely with excitation
    if (x_component == 0.0f && y_component == 0.0f) {
        return x_ekf(0); // Return current EKF estimate if measurement is singular
    }

    return atan2(y_component, x_component);
}

/**
 * @brief Auto-tunes the EXCITATION_PWM_VALUE to center the average phase voltage.
 * Call this function when the motor is stationary.
 */
void autoTuneExcitationPwm() {
    Serial.println("\n--- Starting Auto-Tune for EXCITATION_PWM_VALUE ---");
    Serial.println("  (IMPORTANT: Ensure motor is ABSOLUTELY STATIONARY!)");
    Serial.println("  Monitoring average phase voltage during tuning.");

    // Re-initialize Kalman filter states and covariance matrix for a clean start
    for(int i = 0; i < N_STATES_KF1; ++i) {
        x_hat_kf1(i) = 0.0f;
        for(int j = 0; j < N_STATES_KF1; ++j) {
            P_kf1(i, j) = (i == j) ? 1000.0f : 0.0f; // High initial uncertainty for KF1
        }
    }

    // For EKF during tuning, we can initialize its state to 0 angle, 0 velocity, high uncertainty
    x_ekf(0) = 0.0f; // Start angle at 0 rad
    x_ekf(1) = 0.0f; // Start velocity at 0 rad/s
    for(int i = 0; i < N_STATES_EKF; ++i) {
        for(int j = 0; j < N_STATES_EKF; ++j) {
            P_ekf(i, j) = (i == j) ? 10.0f : 0.0f; // High initial uncertainty for EKF states
        }
    }
    prev_ekf_micros = micros(); // Reset EKF dt timer

    // Reset velocity history
    velocity_history_index = 0;
    velocity_history_filled = false;
    current_rotation_state = STATIONARY;
    consistent_direction_count = 0;

    for (int i = 0; i < PWM_CALIBRATION_ITERATIONS; i++) {
        executeSensingCycle(); // Get current sensor readings with current PWM value (now filtered by KF1)

        int current_avg_adc = (int)((g_valA + g_valB + g_valC) / 3.0f);

        Serial.print("Tune Step "); Serial.print(i + 1);
        Serial.print("\tPWM: "); Serial.print(EXCITATION_PWM_VALUE);
        Serial.print("\tAvg ADC: "); Serial.print(current_avg_adc); // Corrected typo here
        Serial.print("\tPhases (A/B/C): "); Serial.print(g_valA, 2);
        Serial.print("/"); Serial.print(g_valB, 2);
        Serial.print("/"); Serial.println(g_valC, 2);

        if (abs(current_avg_adc - TARGET_QUIESCENT_AVG_ADC) <= PWM_CALIBRATION_TOLERANCE) {
            Serial.println("--- Auto-tune complete! Optimal PWM found. ---");
            break;
        }

        if (current_avg_adc < TARGET_QUIESCENT_AVG_ADC) {
            EXCITATION_PWM_VALUE += PWM_ADJUST_STEP;
        } else {
            EXCITATION_PWM_VALUE -= PWM_ADJUST_STEP;
        }

        if (EXCITATION_PWM_VALUE < 0) EXCITATION_PWM_VALUE = 0;
        if (EXCITATION_PWM_VALUE > 255) EXCITATION_PWM_VALUE = 255;

        delay(50); // Small delay between calibration steps for stability
    }
    Serial.print("Final auto-tuned EXCITATION_PWM_VALUE: "); Serial.println(EXCITATION_PWM_VALUE);
    Serial.println("--------------------------------------------------");
    Serial.println("Resume operation. Rotate motor for data.");
    Serial.println(" "); 
}

void setup() {
    Serial.begin(115200);
    Serial.println("BLDC Encoder Initialized.");

    // Initialize Identity matrices manually for BLA operations
    // For KF1 (3x3 identity)
    I_matrix_kf1 = {
        1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 1.0f
    };
    // For EKF (2x2 identity)
    I_matrix_ekf = {
        1.0f, 0.0f,
        0.0f, 1.0f
    };

    // Initialize KF1 matrices
    F_kf1 = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity
    H_kf1 = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity
    Q_kf1 = {Q_DIAGONAL_VALUE_KF1, 0, 0, 0, Q_DIAGONAL_VALUE_KF1, 0, 0, 0, Q_DIAGONAL_VALUE_KF1};
    R_kf1 = {R_DIAGONAL_VALUE_KF1, 0, 0, 0, R_DIAGONAL_VALUE_KF1, 0, 0, 0, R_DIAGONAL_VALUE_KF1};
    // Initial P_kf1 will be set in autoTuneExcitationPwm()

    // Initialize velocity history
    for (int i = 0; i < VELOCITY_HISTORY_SIZE; i++) {
        velocity_history[i] = 0.0f;
    }

    // ADC setup (faster conversion, internal 1.1V reference)
    analogReference(INTERNAL);
    ADCSRA &= ~((1 << ADPS0) | (1 << ADPS1) | (1 << ADPS2)); // Clear prescaler bits
    ADCSRA |= (1 << ADPS2) | (1 << ADPS0); // Set prescaler to 32 (16MHz / 32 = 500kHz ADC clock)
    ADCSRA |= (1 << ADEN); // Enable ADC

    // INCREASE PWM FREQUENCY FOR PIN 9 (Timer1)
    TCCR1B = (TCCR1B & 0xF8) | 0x01; // Set CS12, CS11, CS10 to 001 for prescaler of 1

    // Common pin setup (default to output LOW, strong sink)
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);

#ifndef DIAGNOSTIC_MODE
    autoTuneExcitationPwm();

    // After tuning, initialize x_hat_kf1 to the target average values
    x_hat_kf1(0) = (float)TARGET_QUIESCENT_AVG_ADC;
    x_hat_kf1(1) = (float)TARGET_QUIESCENT_AVG_ADC;
    x_hat_kf1(2) = (float)TARGET_QUIESCENT_AVG_ADC;
    // Lower uncertainty for KF1 after tuning
    for(int i = 0; i < N_STATES_KF1; ++i) {
        for(int j = 0; j < N_STATES_KF1; ++j) {
            P_kf1(i, j) = (i == j) ? 100.0f : 0.0f;
        }
    }

    // Initialize EKF states based on a first measurement (approximate)
    executeSensingCycle(); // Get initial filtered values
    x_ekf(0) = getInferredElectricalAngle(); // Initial angle guess from current readings
    x_ekf(1) = 0.0f;                         // Initial velocity guess (stationary)

    // Initial uncertainty for EKF (P_ekf)
    P_ekf = {
        1.0f, 0,  // Angle uncertainty (rad^2)
        0, 1.0f   // Velocity uncertainty ((rad/s)^2)
    };
    prev_ekf_micros = micros(); // Start EKF dt timer
#endif
}

void loop() {
#ifdef DIAGNOSTIC_MODE

    executeSensingCycle(); // Get latest filtered phase voltages (g_valA/B/C)

    float inferred_angle = getInferredElectricalAngle(); // Measurement for EKF
    applyEKF(inferred_angle); // Apply EKF to estimate angle and velocity

    // Convert estimated electrical angle to 'Hall steps' for easier comparison
    // Each electrical cycle (2*PI rad) corresponds to 6 Hall steps
    estimated_mechanical_steps = round(x_ekf(0) / (2.0f * PI) * 6.0f);

    // Get rotation direction using improved spurious rotation prevention
    int rotation_direction = getRotationDirection();

    // Get additional diagnostic information
    float avg_velocity = getAverageVelocity();
    float velocity_std = getVelocityStandardDeviation();

    // Output for Serial Plotter:
    // Filtered A, Filtered B, Filtered C, Inferred Angle (deg), EKF Angle (deg), 
    // Raw EKF Velocity (rad/s), Averaged Velocity (rad/s), Velocity Std, 
    // Estimated Mechanical Steps, Rotation Direction
    Serial.print(g_valA, 2);
    Serial.print("\t");
    Serial.print(g_valB, 2);
    Serial.print("\t");
    Serial.print(g_valC, 2);
    Serial.print("\t");
    Serial.print(inferred_angle * 180.0 / PI, 2); // Inferred angle in degrees
    Serial.print("\t");
    Serial.print(x_ekf(0) * 180.0 / PI, 2);       // EKF estimated angle in degrees
    Serial.print("\t");
    Serial.print(x_ 
