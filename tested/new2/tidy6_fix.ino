#include <Arduino.h>
#include <BasicLinearAlgebra.h>

// --- Configuration Namespace ---
namespace Config {
    // Pin Definitions
    constexpr int PHASE_A_PIN = A0;
    constexpr int PHASE_B_PIN = A1;
    constexpr int PHASE_C_PIN = A2;
    constexpr int COMMON_PIN = 9;

    // Motor & Angle Parameters
    constexpr float ELECTRICAL_STEPS_PER_REV = 6.0f; // For mechanical steps calculation
    constexpr float TWO_PI_F = 2.0f * PI;
    constexpr float SQRT_3 = 1.7320508f;

    // ADC & Excitation Parameters
    int EXCITATION_PWM_VALUE = 1; // Calibrated during auto-tuning, mutable
    constexpr int FLOATING_SETTLING_TIME_US = 10;     // Time to settle after excitation/drive

    // New Excitation Parameters for Velocity-Based Pulses & Oversampling
    constexpr int MIN_EXCITATION_PULSE_WIDTH_US = 50;  // Min pulse width (at zero/low velocity)
    constexpr int MAX_EXCITATION_PULSE_WIDTH_US = 1000; // Max pulse width (at high velocity)
    constexpr float MAX_EXPECTED_VELOCITY = 15.0f; // Max velocity in rad/s, for scaling (tune this)

    constexpr int MIN_OVERSAMPLING_SAMPLES = 4;   // Min samples per phase (at high velocity)
    constexpr int MAX_OVERSAMPLING_SAMPLES = 64; // Max samples per phase (at zero/low velocity)

    // Autotuning Parameters
    constexpr int TARGET_QUIESCENT_AVG_ADC = 80;
    constexpr int PWM_ADJUST_STEP = 1;
    constexpr int PWM_CALIBRATION_ITERATIONS = 2000;
    constexpr int PWM_CALIBRATION_TOLERANCE = 1;
    constexpr int AUTOTUNE_PULSE_LENGTH_US = 500; // Fixed pulse length for auto-tuning, mid-range.
    constexpr int AUTOTUNE_OVERSAMPLING = 32;       // Fixed oversampling for auto-tuning.

    // EKF & KF Noise Parameters (Adjusted for potentially better performance)
    constexpr float Q_DIAGONAL_VALUE_KF1 = 0.5f;     // Process noise for KF1 (lower for less trust in process)
    constexpr float R_DIAGONAL_VALUE_KF1 = 200.0f; // Measurement noise for KF1 (higher for more noise in ADC)
    constexpr float KF1_INVERSION_EPSILON = 1e-6f; // Small value to prevent singularity in KF1 S matrix inversion

    constexpr float Q_ANGLE_EKF = 0.0005f; // Process noise for angle (lower for smoother angle)
    constexpr float Q_VEL_EKF = 0.005f;      // Process noise for velocity (lower for smoother velocity)
    constexpr float R_ANGLE_EKF = 0.2f;      // Measurement noise for EKF (higher for more noise in inferred angle)
    constexpr float EKF_INVERSION_EPSILON = 1e-6f; // Small value to prevent singularity in EKF S matrix inversion
    constexpr float EKF_OUTLIER_THRESHOLD = PI; // Max innovation (angle diff) allowed before considering it an outlier

    // Velocity & Locking Logic Parameters
    constexpr int VELOCITY_HISTORY_SIZE = 15; // Increased history for better smoothing/trend detection
    constexpr int ANGLE_HISTORY_SIZE = 7;       // Increased history for better consistency checks
    constexpr float ANGLE_UPDATE_DEADBAND = 0.02f; // Deadband for angle updates (not directly used, consider for future use)
    constexpr float VELOCITY_STEP_THRESHOLD = 0.04f; // Reduced threshold for "is_moving"
    constexpr float MIN_VELOCITY_CONSISTENCY = 0.05f; // Minimum velocity for consistent direction (increased slightly)
    constexpr float VELOCITY_NOISE_THRESHOLD = 0.015f; // Threshold for velocity noise rejection (lower for stricter consistency)
    constexpr unsigned long ZERO_VEL_LOCK_TIMEOUT_MS = 250UL; // Timeout for zero velocity unlock (ms)
    constexpr unsigned long MIN_STABLE_TIME_US = 50000UL; // Minimum time for stable velocity (currently unused but good to keep)

    // Revolution Counting
    constexpr float REVOLUTION_THRESHOLD_FACTOR = 0.75f; // Angle difference factor for revolution counting (PI * factor)
    constexpr int CONSISTENT_DIR_COUNT_THRESHOLD = 5; // How many consistent velocity samples before locking
    constexpr float VELOCITY_TREND_THRESHOLD = 0.008f; // Trend magnitude to consider significant for stopping (adjusted)

    // --- Parameter Validation (Compile-Time Assertions) ---
    static_assert(MIN_EXCITATION_PULSE_WIDTH_US > 0, "MIN_EXCITATION_PULSE_WIDTH_US must be positive.");
    static_assert(MAX_EXCITATION_PULSE_WIDTH_US >= MIN_EXCITATION_PULSE_WIDTH_US, "MAX_EXCITATION_PULSE_WIDTH_US must be >= MIN_EXCITATION_PULSE_WIDTH_US.");
    static_assert(MIN_OVERSAMPLING_SAMPLES > 0, "MIN_OVERSAMPLING_SAMPLES must be positive.");
    static_assert(MAX_OVERSAMPLING_SAMPLES >= MIN_OVERSAMPLING_SAMPLES, "MAX_OVERSAMPLING_SAMPLES must be >= MIN_OVERSAMPLING_SAMPLES.");
    static_assert(MAX_EXPECTED_VELOCITY > 0.0f, "MAX_EXPECTED_VELOCITY must be positive.");
    static_assert(TARGET_QUIESCENT_AVG_ADC >= 0 && TARGET_QUIESCENT_AVG_ADC <= 1023, "TARGET_QUIESCENT_AVG_ADC must be within ADC range.");
    static_assert(PWM_ADJUST_STEP > 0, "PWM_ADJUST_STEP must be positive.");
    static_assert(FLOATING_SETTLING_TIME_US >= 0, "FLOATING_SETTLING_TIME_US must be non-negative.");
    static_assert(KF1_INVERSION_EPSILON >= 0.0f, "KF1_INVERSION_EPSILON must be non-negative.");
    static_assert(EKF_INVERSION_EPSILON >= 0.0f, "EKF_INVERSION_EPSILON must be non-negative.");
    static_assert(EKF_OUTLIER_THRESHOLD > 0.0f, "EKF_OUTLIER_THRESHOLD must be positive.");
}

// --- Global Variables ---
#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;               // KF1 state: [Phase A, Phase B, Phase C] (Offset Corrected)
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1;

#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf;              // EKF state: [Angle, Velocity]
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf;
BLA::Matrix<1, 1> R_ekf;

// Store initial P_ekf for potential reset
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> INITIAL_P_EKF;
float last_valid_inferred_angle = 0.0f; // To store last known good angle

// --- FIX: Add a global variable to store the calibrated quiescent offset ---
int quiescent_offset_adc = Config::TARGET_QUIESCENT_AVG_ADC;

// --- Enhanced State Tracking History ---
struct HistoryEntry {
    float angle;
    float velocity;
    unsigned long timestamp;
};

HistoryEntry velocityHistory[Config::VELOCITY_HISTORY_SIZE];
HistoryEntry angleHistory[Config::ANGLE_HISTORY_SIZE];
int velocityHistoryIndex = 0;
int angleHistoryIndex = 0;
bool velocityHistoryFull = false;
bool angleHistoryFull = false;

// --- EKF and Motor State Variables ---
bool ekfLocked = false;
unsigned long lastNonZeroVelMicros = 0;
unsigned long prev_ekf_micros = 0;
long revolution_count = 0;
long estimated_mechanical_steps = 0;
float unwrapped_electrical_angle = 0.0f;

// --- Additional Tracking and Logic Variables ---
float smoothed_velocity = 0.0f;
float velocity_trend = 0.0f;
unsigned long last_direction_change_micros = 0;
int consistent_direction_count = 0;
float last_stable_velocity = 0.0f; // Stores the last velocity when consistently moving
bool is_stopping = false; // Flag indicating if rotation is actively stopping
bool prev_is_moving = false; // To detect transitions from moving to not moving

// --- Function Prototypes ---
void initializeMatrices();
void autoTuneExcitationPwm();
int getExcitationPulseWidth();     // Dynamic pulse width
int getNumOversamplingSamples();   // Dynamic oversampling
void executeSensingCycle();
void applyMultivariateKalmanFilter(const int readings[3], int numSamples); // numSamples passed to KF
float getInferredElectricalAngle();
float normalizeAngle(float angle);
void applyEKF(float measured_angle);
void logData();
int checkedAnalogRead(int pin);

// New functions for improved tracking
void updateVelocityHistory(float velocity, unsigned long timestamp);
void updateAngleHistory(float angle, unsigned long timestamp);
float getSmoothedVelocity();
float getVelocityTrend();
bool isVelocityConsistent(float velocity_to_check);
bool isRotationStopping();
void updateRevolutionCount(float angle_diff);

// --- Setup Function ---
void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10); // Wait for serial port to connect.

    Serial.println(F("Initializing Enhanced BLDC Sensorless Encoder..."));

    initializeMatrices(); // Initialize BLA matrices
    analogReference(INTERNAL); // Use internal 1.1V reference for ADC for stability
    // Configure ADC for faster conversion (Prescaler to 32 -> 16MHz/32 = 500kHz ADC clock)
    ADCSRA = (ADCSRA & ~0x07) | 0x05; // Set prescaler to 32 (0x05 for 32, 0x04 for 16)
    pinMode(Config::COMMON_PIN, OUTPUT); // Ensure common pin is in OUTPUT mode from the start
    digitalWrite(Config::COMMON_PIN, LOW); // Start with common pin low

    autoTuneExcitationPwm(); // Calibrate excitation PWM

    // --- FIX: Initialize KF1 state to zero, as it now represents the offset-corrected signal (BEMF) ---
    x_hat_kf1.Fill(0.0f);
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) P_kf1(i, i) = 100.0f; // High initial uncertainty

    // Initialize history arrays
    for (int i = 0; i < Config::VELOCITY_HISTORY_SIZE; i++) {
        velocityHistory[i] = {0.0f, 0.0f, 0UL};
    }
    for (int i = 0; i < Config::ANGLE_HISTORY_SIZE; i++) {
        angleHistory[i] = {0.0f, 0.0f, 0UL};
    }

    // Perform initial sensing cycle to get first inferred angle
    executeSensingCycle();
    x_ekf(0) = getInferredElectricalAngle(); // Initial angle from KF1
    last_valid_inferred_angle = x_ekf(0); // Store initial valid angle
    x_ekf(1) = 0.0f;                       // Initial velocity is zero
    P_ekf.Fill(0.0f);
    P_ekf(0, 0) = 1.0f; // Initial angle uncertainty
    P_ekf(1, 1) = 1.0f; // Initial velocity uncertainty
    INITIAL_P_EKF = P_ekf; // Save for potential reset

    prev_ekf_micros = micros();
    unwrapped_electrical_angle = x_ekf(0);

    Serial.println(F("Enhanced system ready. Rotate motor to begin estimation."));
}

// --- Main Loop Function ---
void loop() {
    executeSensingCycle();
    float inferred_angle = getInferredElectricalAngle(); // This function now handles NaN/Inf internally

    // Check if the inferred angle is valid before passing to EKF
    if (isnan(inferred_angle) || isinf(inferred_angle)) {
        Serial.println(F("CRITICAL ERROR: Invalid inferred angle (NaN/Inf) detected. Skipping EKF update."));
        // Potentially reset EKF here or revert to last good state if persistent
        // For now, EKF will continue prediction, but measurement update will effectively be skipped.
    } else {
        last_valid_inferred_angle = inferred_angle; // Update last known good angle
        applyEKF(inferred_angle);
    }


    // Calculate mechanical steps based on unwrapped electrical angle
    estimated_mechanical_steps = round(unwrapped_electrical_angle / Config::TWO_PI_F * Config::ELECTRICAL_STEPS_PER_REV);

    logData(); // Log data to serial for debugging/monitoring

    // Handle serial commands for re-tuning
    if (Serial.available() > 0) {
        char command = toupper(Serial.read());
        while (Serial.available() > 0) Serial.read(); // Clear serial buffer
        if (command == 'T') {
            Serial.println(F("\nRe-tuning requested..."));
            autoTuneExcitationPwm(); // Re-calibrate PWM
            executeSensingCycle();   // Get new initial readings
            x_ekf(0) = getInferredElectricalAngle(); // Reset EKF angle
            last_valid_inferred_angle = x_ekf(0); // Update last good angle
            x_ekf(1) = 0.0f;                       // Reset EKF velocity
            P_ekf = INITIAL_P_EKF; // Reset EKF uncertainty to initial value
            prev_ekf_micros = micros();
            revolution_count = 0;
            unwrapped_electrical_angle = x_ekf(0);
            ekfLocked = false; // Unlock EKF after re-tuning

            // Reset history and state tracking variables
            velocityHistoryIndex = angleHistoryIndex = 0;
            velocityHistoryFull = angleHistoryFull = false;
            smoothed_velocity = velocity_trend = 0.0f;
            consistent_direction_count = 0;
            is_stopping = false;
            prev_is_moving = false;

            Serial.println(F("Re-tuning complete. Enhanced system ready."));
        }
    }
}

// --- Function Definitions ---

void initializeMatrices() {
    // Initialize identity matrices
    I_matrix_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) I_matrix_kf1(i, i) = 1.0f;

    I_matrix_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_EKF; i++) I_matrix_ekf(i, i) = 1.0f;

    // Initialize KF1 process and measurement noise covariances
    Q_kf1.Fill(0.0f);
    R_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) {
        Q_kf1(i, i) = Config::Q_DIAGONAL_VALUE_KF1;
        R_kf1(i, i) = Config::R_DIAGONAL_VALUE_KF1;
    }

    // Initialize EKF process and measurement noise covariances
    Q_ekf.Fill(0.0f);
    Q_ekf(0, 0) = Config::Q_ANGLE_EKF; // Angle noise
    Q_ekf(1, 1) = Config::Q_VEL_EKF;   // Velocity noise

    R_ekf(0, 0) = Config::R_ANGLE_EKF; // Angle measurement noise
}

// Auto-tuning for the base PWM value (for quiescent state)
void autoTuneExcitationPwm() {
    int sum = 0;
    Serial.print(F("Auto-tuning PWM. Target ADC: ")); Serial.println(Config::TARGET_QUIESCENT_AVG_ADC);
    Serial.print(F("Using fixed pulse length: ")); Serial.print(Config::AUTOTUNE_PULSE_LENGTH_US); Serial.println(F(" us"));
    Serial.print(F("Using fixed oversampling: ")); Serial.print(Config::AUTOTUNE_OVERSAMPLING); Serial.println(F(" samples"));

    for (int iter = 0; iter < Config::PWM_CALIBRATION_ITERATIONS; iter++) {
        sum = 0;
        for (int i = 0; i < Config::AUTOTUNE_OVERSAMPLING; i++) {
            digitalWrite(Config::COMMON_PIN, LOW); // Ensure low before excitation
            analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE); // Excite common pin
            delayMicroseconds(Config::AUTOTUNE_PULSE_LENGTH_US); // Fixed pulse for auto-tuning

            // --- REVISED SENSING LOGIC ---
            digitalWrite(Config::COMMON_PIN, LOW); // Actively drive common LOW for sensing
            // Pin remains in OUTPUT mode. This provides a ground reference for BEMF.
            // --- END REVISED SENSING LOGIC ---
            
            delayMicroseconds(Config::FLOATING_SETTLING_TIME_US); // Allow voltages to settle

            int reading = checkedAnalogRead(Config::PHASE_A_PIN); // Measure one phase
            sum += reading;
        }
        int avg = sum / Config::AUTOTUNE_OVERSAMPLING;

        if (abs(avg - Config::TARGET_QUIESCENT_AVG_ADC) <= Config::PWM_CALIBRATION_TOLERANCE) {
            Serial.print(F("PWM calibration converged in ")); Serial.print(iter + 1); Serial.println(F(" iterations."));
            break; // Converged
        }

        if (avg < Config::TARGET_QUIESCENT_AVG_ADC) {
            Config::EXCITATION_PWM_VALUE += Config::PWM_ADJUST_STEP;
        } else {
            Config::EXCITATION_PWM_VALUE -= Config::PWM_ADJUST_STEP;
        }
        Config::EXCITATION_PWM_VALUE = constrain(Config::EXCITATION_PWM_VALUE, 0, 255);
        if (Config::EXCITATION_PWM_VALUE == 0 && avg < Config::TARGET_QUIESCENT_AVG_ADC) {
             Serial.println(F("Warning: PWM hit 0, cannot reach target ADC. Adjust TARGET_QUIESCENT_AVG_ADC or hardware."));
             break; // Cannot go lower, break to prevent infinite loop
        }
        if (Config::EXCITATION_PWM_VALUE == 255 && avg > Config::TARGET_QUIESCENT_AVG_ADC) {
             Serial.println(F("Warning: PWM hit 255, cannot reach target ADC. Adjust TARGET_QUIESCENT_AVG_ADC or hardware."));
             break; // Cannot go higher, break to prevent infinite loop
        }
    }
    Serial.print(F("Calibrated base PWM value: "));
    Serial.println(Config::EXCITATION_PWM_VALUE);

    // --- FIX: Set the global quiescent offset to the target value after tuning ---
    quiescent_offset_adc = Config::TARGET_QUIESCENT_AVG_ADC;
    Serial.print(F("Quiescent ADC offset set to: "));
    Serial.println(quiescent_offset_adc);

    digitalWrite(Config::COMMON_PIN, LOW); // Ensure pin is low after calibration
}

// Determine excitation pulse width based on velocity
int getExcitationPulseWidth() {
    float current_velocity = fabsf(x_ekf(1)); // Absolute velocity from EKF

    // Normalize velocity between 0 and 1, where 0 is zero velocity, 1 is MAX_EXPECTED_VELOCITY
    float normalized_velocity = constrain(current_velocity / Config::MAX_EXPECTED_VELOCITY, 0.0f, 1.0f);

    // Interpolate pulse width:
    // At normalized_velocity = 0 (low speed), pulse_width = MIN_EXCITATION_PULSE_WIDTH_US
    // At normalized_velocity = 1 (high speed), pulse_width = MAX_EXCITATION_PULSE_WIDTH_US
    // This makes pulse width directly proportional to velocity.
    int pulse_width = map(normalized_velocity * 1000, 0, 1000,
                          Config::MIN_EXCITATION_PULSE_WIDTH_US, Config::MAX_EXCITATION_PULSE_WIDTH_US);

    return constrain(pulse_width, Config::MIN_EXCITATION_PULSE_WIDTH_US, Config::MAX_EXCITATION_PULSE_WIDTH_US);
}

// Determine number of oversampling samples based on velocity
int getNumOversamplingSamples() {
    float current_velocity = fabsf(x_ekf(1)); // Absolute velocity from EKF

    // Normalize velocity between 0 and 1
    float normalized_velocity = constrain(current_velocity / Config::MAX_EXPECTED_VELOCITY, 0.0f, 1.0f);

    // Interpolate oversampling samples:
    // At normalized_velocity = 0 (low speed), samples = MAX_OVERSAMPLING_SAMPLES
    // At normalized_velocity = 1 (high speed), samples = MIN_OVERSAMPLING_SAMPLES
    // This makes oversampling inversely proportional to velocity.
    int num_samples = map(normalized_velocity * 1000, 0, 1000,
                          Config::MAX_OVERSAMPLING_SAMPLES, Config::MIN_OVERSAMPLING_SAMPLES);

    return constrain(num_samples, Config::MIN_OVERSAMPLING_SAMPLES, Config::MAX_OVERSAMPLING_SAMPLES);
}


void executeSensingCycle() {
    int readings[3] = {0}; // Raw ADC sums for each phase
    int phasePins[] = {Config::PHASE_A_PIN, Config::PHASE_B_PIN, Config::PHASE_C_PIN};
    
    // Determine dynamic pulse width and oversampling
    int current_excitation_pulse_width_us = getExcitationPulseWidth();
    int current_num_oversampling_samples = getNumOversamplingSamples();

    // Oversampling is applied here to reduce random measurement noise before Kalman filtering.
    // The Kalman filter then uses these averaged measurements, along with its dynamic model,
    // to estimate the true state, leveraging inter-phase relationships and system physics.
    for (int i = 0; i < current_num_oversampling_samples; i++) {
        digitalWrite(Config::COMMON_PIN, LOW); // Ensure low before excitation
        analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE); // Excite common pin
        delayMicroseconds(current_excitation_pulse_width_us); // Apply excitation pulse

        // --- REVISED SENSING LOGIC ---
        // After excitation, actively drive the common pin LOW (to ground).
        // This provides a stable reference point for the BEMF measurement on the other phases.
        digitalWrite(Config::COMMON_PIN, LOW);
        // Pin remains in OUTPUT mode.
        // --- END REVISED SENSING LOGIC ---

        delayMicroseconds(Config::FLOATING_SETTLING_TIME_US); // Allow voltages to settle relative to grounded common

        // Read all three phases while common is low
        readings[0] += checkedAnalogRead(phasePins[0]);
        readings[1] += checkedAnalogRead(phasePins[1]);
        readings[2] += checkedAnalogRead(phasePins[2]);
    }

    // Average the readings
    for (int i = 0; i < 3; i++) {
        readings[i] /= current_num_oversampling_samples;
        // --- FIX: Subtract the calibrated offset to normalize the reading around zero ---
        readings[i] -= quiescent_offset_adc;
    }

    applyMultivariateKalmanFilter(readings, current_num_oversampling_samples);
}

void applyMultivariateKalmanFilter(const int readings[3], int numSamples) {
    BLA::Matrix<N_STATES_KF1, 1> z; // Measurement vector
    for (int i = 0; i < 3; i++) z(i) = readings[i];

    // Kalman Filter equations
    BLA::Matrix<N_STATES_KF1, 1> y = z - x_hat_kf1; // Innovation
    
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S = P_kf1 + R_kf1; // Innovation covariance
    
    // Add small epsilon to diagonal to prevent singularity for S matrix inversion (KF1)
    for (int i = 0; i < N_STATES_KF1; ++i) {
        S(i,i) += Config::KF1_INVERSION_EPSILON;
    }

    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> K = P_kf1 * BLA::Inverse(S); // Kalman gain

    x_hat_kf1 += K * y; // Update state estimate
    P_kf1 = (I_matrix_kf1 - K) * P_kf1 + Q_kf1; // Update error covariance
}

// ... (The rest of the file remains unchanged as the fix is upstream) ...
