#include <Arduino.h>
#include <BasicLinearAlgebra.h>

namespace Config {
    constexpr int PHASE_A_PIN = A0, PHASE_B_PIN = A1, PHASE_C_PIN = A2, COMMON_PIN = 9;
    constexpr float ELECTRICAL_STEPS_PER_REV = 1.0f;
    
    // Define PI and TWO_PI_F for clarity and consistency
    const float PI_F = 3.14159265358979323846f;
    constexpr float TWO_PI_F = 2.0f * PI_F; // PI_F is now used
    constexpr float SQRT_3 = 1.7320508f;
    
    // Original NUM_SAMPLES is now a default/initial value
    constexpr int DEFAULT_NUM_SAMPLES = 16; 

    // New parameters for dynamic oversampling
    constexpr int MIN_NUM_SAMPLES = 16;  // Minimum number of samples
    constexpr int MAX_NUM_SAMPLES = 16;  // Maximum number of samples
    // **MODIFIED THRESHOLDS FOR VELOCITY**
    constexpr float LOW_VELOCITY_THRESHOLD = 10.0f; // Velocity below which more samples are taken (rad/s)
    constexpr float HIGH_VELOCITY_THRESHOLD = 10.0f; // Velocity above which fewer samples are taken (rad/s)

    constexpr float ANGLE_UPDATE_DEADBAND = 0.0002f, VELOCITY_STEP_THRESHOLD = 0.0001f; // Reduced threshold
    constexpr unsigned long ZERO_VEL_LOCK_TIMEOUT = 200000UL; // Reduced timeout
    int EXCITATION_PWM_VALUE = 1;
    constexpr int EXCITATION_PULSE_WIDTH_US = 32, FLOATING_SETTLING_TIME_US = 1;
    constexpr int TARGET_QUIESCENT_AVG_ADC = 128, PWM_ADJUST_STEP = 1, PWM_CALIBRATION_ITERATIONS = 2000, PWM_CALIBRATION_TOLERENCE = 1;
    constexpr float Q_DIAGONAL_VALUE_KF1 = 5.0f, R_DIAGONAL_VALUE_KF1 = 80.0f;
    constexpr float Q_ANGLE_EKF = 1.1f, Q_VEL_EKF = 80.1f; 
    constexpr float R_ANGLE_EKF = 1.1f; // Base R_ekf value
    constexpr float R_ANGLE_EKF_MAX = 10.0f; // Max R_ekf value when phase consistency is low

    // New parameters for improved tracking
    constexpr int VELOCITY_HISTORY_SIZE = 4;
    constexpr int ANGLE_HISTORY_SIZE = 4;
    constexpr int PHASE_HISTORY_SIZE = 8; // Increased history for KF filtered phases for better regression
    constexpr float MIN_VELOCITY_CONSISTENCY = 0.0003f; // Minimum velocity for consistent direction
    constexpr float VELOCITY_NOISE_THRESHOLD = 0.006f; // Threshold for velocity noise rejection
    constexpr float ANGLE_HISTORY_NOISE_THRESHOLD = 0.0005f; // Threshold for angle consistency from phase history (radians) - now less critical with R-squared
    constexpr unsigned long MIN_STABLE_TIME_US = 500UL; // Minimum time for stable velocity
    constexpr float PHASE_CONSISTENCY_R_SQUARED_THRESHOLD = 0.85f; // R-squared threshold for considering phase history consistent
}

#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1; // KF state estimate (filtered ADC readings)
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1; // KF covariance matrices

#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf; // EKF state estimate (angle, velocity)
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf; // EKF covariance matrices
BLA::Matrix<1, 1> R_ekf; // EKF measurement noise covariance

// Enhanced state tracking structs
struct HistoryEntry {
    float angle;
    float velocity;
    unsigned long timestamp;
};

// New struct for phase history
struct PhaseHistoryEntry {
    float phaseA, phaseB, phaseC;
    unsigned long timestamp;
};

// History buffers
HistoryEntry velocityHistory[Config::VELOCITY_HISTORY_SIZE];
HistoryEntry angleHistory[Config::ANGLE_HISTORY_SIZE];
PhaseHistoryEntry phaseHistory[Config::PHASE_HISTORY_SIZE]; // New phase history buffer
int velocityHistoryIndex = 0, angleHistoryIndex = 0, phaseHistoryIndex = 0; 
bool velocityHistoryFull = false, angleHistoryFull = false, phaseHistoryFull = false;

// Global variables for the new rotation tracking system
float previous_raw_electrical_angle = 0.0f; // Stores the previous raw angle from KF for unwrapping
bool first_angle_reading = true; // Flag to initialize unwrapped angle on first read

// Main angle tracking variables
float unwrapped_electrical_angle = 0.0f; // This is the continuous, unbounded angle (EKF's angle, unwrapped)
long revolution_count = 0;
long estimated_mechanical_steps = 0;

// EKF and consistency flags
bool ekfLocked = false; // EKF lock indicates confidence
unsigned long lastNonZeroVelMicros = 0, prev_ekf_micros = 0;

float smoothed_velocity = 0.0f;
float velocity_trend = 0.0f;
unsigned long last_direction_change = 0; // Not used in current version, but kept
int consistent_direction_count = 0; // Not fully utilized in current version, but kept for future use
float last_stable_velocity = 0.0f; // Not used in current version, but kept
bool is_stopping = false;
bool is_consistent_phases = false; // Flag to indicate phase consistency based on R-squared

// Variables for the "Novel Method" candidate simulation
float historical_inferred_velocity = 0.0f; // Velocity from linear regression on phase history
float historical_inferred_angle = 0.0f;    // Angle projected from linear regression on phase history
float phase_history_consistency_score = 0.0f; // R-squared of the linear regression

// Function declarations
void initializeMatrices();
void autoTuneExcitationPwm();
void executeSensingCycle();
float getInferredElectricalAngle();
float normalizeAngle(float);
void applyEKF(float);
void logData();
int checkedAnalogRead(int);
void applyMultivariateKalmanFilter(const long readings[3]); 
void updateVelocityHistory(float velocity, unsigned long timestamp);
void updateAngleHistory(float angle, unsigned long timestamp);
void updatePhaseHistory(float phaseA, float phaseB, float phaseC, unsigned long timestamp);

// Struct and function for the novel historical phase inference
struct HistoricalPhaseStats {
    float inferredVelocity;
    float inferredAngleAtCurrentTime; // Angle projected to the current time
    float consistencyScore;           // e.g., R-squared of the linear regression
    bool isValid;                     // Indicates if enough data for valid stats
};
HistoricalPhaseStats getHistoricalPhaseStats(); // New function for phase history analysis

bool isVelocityConsistent();
bool isRotationStopping(); // Determines if rotation is stopping based on velocity trend

void updateUnwrappedElectricalAngle(float current_raw_electrical_angle);

int getDynamicNumSamples(float current_velocity);

// Kahan Summation for precision in floating-point sums
void kahanSum(float newValue, float &sum, float &compensation) {
    float y = newValue - compensation;
    float t = sum + y;
    compensation = (t - sum) - y; // Calculate the error (lost low-order bits)
    sum = t; 
}

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println(F("Initializing Enhanced BLDC Sensorless Encoder..."));
    initializeMatrices();
    analogReference(INTERNAL);
    // Configure ADC prescaler for faster readings (e.g., 32 for 16MHz clock -> 38.5kHz)
    ADCSRA = (ADCSRA & ~0x07) | 0x05; // Set prescaler to 32 (16MHz/32 = 500kHz, 13 cycles/read = 38.5k samples/sec)
    // Configure Timer1 for 1us resolution for delayMicroseconds (assuming 16MHz Arduino, TCCR1B = 0x01 means prescaler 1)
    // This also sets Timer1's prescaler for PWM on pins 9 and 10 to 1, resulting in 62.5 kHz PWM frequency.
    TCCR1B = (TCCR1B & 0b11111000) | 0x01; // Set Timer1 prescaler to 1 (No prescaling for 16MHz)

    pinMode(Config::COMMON_PIN, OUTPUT); // Ensure COMMON_PIN is set as OUTPUT initially
    digitalWrite(Config::COMMON_PIN, LOW); // Start with common pin low
    autoTuneExcitationPwm();

    x_hat_kf1.Fill((float)Config::TARGET_QUIESCENT_AVG_ADC);
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) P_kf1(i, i) = 100.0f;

    // Initialize history buffers
    for (int i = 0; i < Config::VELOCITY_HISTORY_SIZE; i++) {
        velocityHistory[i] = {0.0f, 0.0f, 0UL};
    }
    for (int i = 0; i < Config::ANGLE_HISTORY_SIZE; i++) {
        angleHistory[i] = {0.0f, 0.0f, 0UL};
    }
    for (int i = 0; i < Config::PHASE_HISTORY_SIZE; i++) {
        phaseHistory[i] = {0.0f, 0.0f, 0.0f, 0UL};
    }

    Serial.println(F("Enhanced system ready. Rotate motor to begin estimation."));
}

void loop() {
    // 1. Execute sensing cycle and update KF states
    executeSensingCycle(); 

    // 2. Get the instantaneous electrical angle from KF filtered phases
    float current_raw_angle_from_kf = getInferredElectricalAngle();

    // 3. Update the continuous, unwrapped electrical angle using the new KF angle
    updateUnwrappedElectricalAngle(current_raw_angle_from_kf);

    // 4. Apply EKF: EKF now filters the current wrapped angle and estimates velocity
    applyEKF(current_raw_angle_from_kf); 

    // 5. Calculate mechanical steps from the unwrapped angle
    // Each electrical revolution (2*PI) corresponds to ELECTRICAL_STEPS_PER_REV mechanical steps.
    estimated_mechanical_steps = round(unwrapped_electrical_angle / Config::TWO_PI_F * Config::ELECTRICAL_STEPS_PER_REV);
    revolution_count = round(unwrapped_electrical_angle / Config::TWO_PI_F); // Full electrical revolutions

    // 6. Log data for Serial Plotter/Monitor
    logData();

    // Handle serial commands for re-tuning
    if (Serial.available() > 0) {
        char command = toupper(Serial.read());
        while (Serial.available() > 0) Serial.read(); // Clear buffer
        if (command == 'T') {
            Serial.println(F("\nRe-tuning requested..."));
            autoTuneExcitationPwm();
            
            // Re-perform initial sensing after auto-tuning
            int initial_num_samples = Config::DEFAULT_NUM_SAMPLES; 
            long readings_initial[3] = {0};
            for (int i = 0; i < initial_num_samples; i++) {
                pinMode(Config::COMMON_PIN, OUTPUT);
                digitalWrite(Config::COMMON_PIN, LOW);
                TCNT1 = 0;
                analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
                delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);
                analogWrite(Config::COMMON_PIN, 0);
                // The original code had pinMode(Config::COMMON_PIN, INPUT); delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
                // These were commented out in the user provided code and replaced implicitly with just analogRead
                // For consistency with sensing cycle, removed them from here as well.
                readings_initial[0] += checkedAnalogRead(Config::PHASE_A_PIN);
                readings_initial[1] += checkedAnalogRead(Config::PHASE_B_PIN);
                readings_initial[2] += checkedAnalogRead(Config::PHASE_C_PIN);
            }
            for (int i = 0; i < 3; i++) readings_initial[i] /= initial_num_samples;
            applyMultivariateKalmanFilter(readings_initial);

            float new_initial_raw_angle = getInferredElectricalAngle();
            x_ekf(0) = new_initial_raw_angle;
            x_ekf(1) = 0.0f;
            P_ekf(0, 0) = 1.0f;
            P_ekf(1, 1) = 1.0f;
            prev_ekf_micros = micros();

            // Reset unwrapped angle tracking
            first_angle_reading = true; // Re-initialize unwrapped angle on next update
            updateUnwrappedElectricalAngle(new_initial_raw_angle); // Re-sets unwrapped_electrical_angle and previous_raw_electrical_angle
            
            revolution_count = 0;
            estimated_mechanical_steps = 0; // Reset step counter

            // Reset history buffers
            velocityHistoryIndex = angleHistoryIndex = phaseHistoryIndex = 0;
            velocityHistoryFull = angleHistoryFull = phaseHistoryFull = false;
            
            // Reset EKF confidence flags
            smoothed_velocity = velocity_trend = 0.0f;
            consistent_direction_count = 0;
            is_stopping = false;
            is_consistent_phases = false; // Reset phase consistency flag
            ekfLocked = false; // Reset lock state

            Serial.println(F("Re-tuning complete. Enhanced system ready."));
        }
    }
}

void initializeMatrices() {
    I_matrix_kf1.Fill(0.0f);
    I_matrix_ekf.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) I_matrix_kf1(i, i) = 1.0f;
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
    // R_ekf initial value will be set dynamically in applyEKF
    R_ekf(0, 0) = Config::R_ANGLE_EKF; 
}

void autoTuneExcitationPwm() {
    Serial.println(F("\n--- Starting Auto-Tune ---"));
    Serial.println(F("Ensure motor is stationary!"));

    // Reset filters to a high-uncertainty state
    x_hat_kf1.Fill(0.0f);
    P_kf1.Fill(0.0f);
    for (int i = 0; i < N_STATES_KF1; i++) P_kf1(i, i) = 1000.0f;

    for (int iteration = 0; iteration < Config::PWM_CALIBRATION_ITERATIONS; iteration++) {
        // Use a fixed number of samples for auto-tuning to ensure consistent measurement
        int auto_tune_samples = Config::DEFAULT_NUM_SAMPLES; // Or a specific count like 32/64 for max accuracy during tune
        long readings_tune[3] = {0};
        for (int i = 0; i < auto_tune_samples; i++) {
            pinMode(Config::COMMON_PIN, OUTPUT);
            digitalWrite(Config::COMMON_PIN, LOW);
            TCNT1 = 0;
            analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
            delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);
            analogWrite(Config::COMMON_PIN, 0);
            readings_tune[0] += checkedAnalogRead(Config::PHASE_A_PIN);
            readings_tune[1] += checkedAnalogRead(Config::PHASE_B_PIN);
            readings_tune[2] += checkedAnalogRead(Config::PHASE_C_PIN);
        }
        for (int i = 0; i < 3; i++) readings_tune[i] /= auto_tune_samples;
        applyMultivariateKalmanFilter(readings_tune);
        
        int current_avg_adc = static_cast<int>((x_hat_kf1(0) + x_hat_kf1(1) + x_hat_kf1(2)) / 3.0f);
        int error = current_avg_adc - Config::TARGET_QUIESCENT_AVG_ADC;

        if (abs(error) <= Config::PWM_CALIBRATION_TOLERENCE) {
            Serial.println(F("--- Auto-tune complete! ---"));
            break;
        }

        Config::EXCITATION_PWM_VALUE += (error < 0) ? Config::PWM_ADJUST_STEP : -Config::PWM_ADJUST_STEP;
        Config::EXCITATION_PWM_VALUE = constrain(Config::EXCITATION_PWM_VALUE, 1, 255); // PWM must be > 0

        delay(5); // Short delay for stability
    }
    Serial.print(F("Final Tuned PWM Value: ")); Serial.println(Config::EXCITATION_PWM_VALUE);
}

// Function for dynamic oversampling - now takes velocity
int getDynamicNumSamples(float current_velocity) {
    float abs_velocity = fabsf(current_velocity); // Use absolute velocity
    
    if (abs_velocity < Config::LOW_VELOCITY_THRESHOLD) {
        return Config::MAX_NUM_SAMPLES; // High oversampling for very low velocities (near stop)
    } else if (abs_velocity > Config::HIGH_VELOCITY_THRESHOLD) {
        return Config::MIN_NUM_SAMPLES; // Low oversampling for high velocities
    } else {
        // Linearly interpolate between MIN_NUM_SAMPLES and MAX_NUM_SAMPLES
        // based on where abs_velocity falls between LOW_VELOCITY_THRESHOLD and HIGH_VELOCITY_THRESHOLD.
        // As velocity increases, num_samples decreases.
        float normalized_velocity = (abs_velocity - Config::LOW_VELOCITY_THRESHOLD) /
                                    (Config::HIGH_VELOCITY_THRESHOLD - Config::LOW_VELOCITY_THRESHOLD);
        
        // Invert the interpolation since higher velocity means fewer samples
        int num_samples = Config::MAX_NUM_SAMPLES -
                          static_cast<int>(normalized_velocity * (Config::MAX_NUM_SAMPLES - Config::MIN_NUM_SAMPLES));
        return constrain(num_samples, Config::MIN_NUM_SAMPLES, Config::MAX_NUM_SAMPLES);
    }
}

void executeSensingCycle() {
    long readings[3] = {0};
    
    // Use the EKF's estimated velocity for dynamic oversampling
    float current_velocity = fabsf(x_ekf(1)); // Get the absolute velocity from the EKF state
    
    int num_samples_to_take = getDynamicNumSamples(current_velocity);

    for (int i = 0; i < num_samples_to_take; i++) {
        // Step 1: Ensure COMMON_PIN is in OUTPUT mode and low before starting PWM
        pinMode(Config::COMMON_PIN, OUTPUT);
        digitalWrite(Config::COMMON_PIN, LOW); 

        // Step 2: Reset Timer1 counter (TCNT1) for precise synchronization.
        TCNT1 = 0; 

        // Step 3: Apply PWM excitation signal
        analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);

        // Step 4: Wait for the specified excitation pulse width.
        delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);

        // Step 5: Stop PWM. Common pin remains OUTPUT LOW (not floating as in original provided code)
        analogWrite(Config::COMMON_PIN, 0); 
        // The original user code commented out pinMode(Config::COMMON_PIN, INPUT); and delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
        // This implies the readings are taken while the common pin is outputting 0V after excitation pulse.

        // Step 6: Perform ADC readings while phases are in the measurement state (floating or pulled low after excitation).
        readings[0] += checkedAnalogRead(Config::PHASE_A_PIN);
        readings[1] += checkedAnalogRead(Config::PHASE_B_PIN);
        readings[2] += checkedAnalogRead(Config::PHASE_C_PIN);
    }
    // Average the readings over the number of samples
    for (int i = 0; i < 3; i++) readings[i] /= num_samples_to_take;
    // Apply the Kalman Filter to the averaged readings
    applyMultivariateKalmanFilter(readings);
    
    // Ensure common pin is low after the cycle for the rest of the time
    pinMode(Config::COMMON_PIN, OUTPUT);
    digitalWrite(Config::COMMON_PIN, LOW); 
}

void applyMultivariateKalmanFilter(const long readings[3]) {
    BLA::Matrix<N_STATES_KF1, 1> z, y;
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S; // S must be N_STATES_KF1 x N_STATES_KF1 (3x3)

    for (int i = 0; i < 3; i++) z(i) = readings[i];
    y = z - x_hat_kf1; // Measurement residual
    S = P_kf1 + R_kf1; // Innovation covariance
    auto K = P_kf1 * BLA::Inverse(S); // Kalman Gain (Inverse will now work on 3x3 S)
    x_hat_kf1 += K * y; // Updated state estimate
    P_kf1 = (I_matrix_kf1 - K) * P_kf1 + Q_kf1; // Updated covariance

    // Update KF phase history
    updatePhaseHistory(x_hat_kf1(0), x_hat_kf1(1), x_hat_kf1(2), micros());
}

// Calculates electrical angle from KF filtered phase readings using Clarke transform
float getInferredElectricalAngle() {
    float a = x_hat_kf1(0), b = x_hat_kf1(1), c = x_hat_kf1(2);
    // Subtract ADC offset before Clarke Transform for better accuracy
    a -= Config::TARGET_QUIESCENT_AVG_ADC;
    b -= Config::TARGET_QUIESCENT_AVG_ADC;
    c -= Config::TARGET_QUIESCENT_AVG_ADC;
    
    // Alpha-Beta transform (Clarke Transform)
    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta = (2.0f / 3.0f) * ((b - c) / Config::SQRT_3);
    return atan2f(beta, alpha); // Returns angle in range (-PI, PI]
}

// Function for continuous angle unwrapping
void updateUnwrappedElectricalAngle(float current_raw_electrical_angle) {
    if (first_angle_reading) {
        unwrapped_electrical_angle = current_raw_electrical_angle;
        previous_raw_electrical_angle = current_raw_electrical_angle;
        first_angle_reading = false;
        return;
    }

    float angle_diff = current_raw_electrical_angle - previous_raw_electrical_angle;

    // Handle angle wrap-around directly
    // If the difference is greater than PI (e.g., crossing from +PI to -PI),
    // it means we've wrapped around counter-clockwise. Add 2*PI to make it positive.
    if (angle_diff > Config::PI_F) {
        angle_diff -= Config::TWO_PI_F;
    } 
    // If the difference is less than -PI (e.g., crossing from -PI to +PI),
    // it means we've wrapped around clockwise. Subtract 2*PI to make it negative.
    else if (angle_diff < -Config::PI_F) {
        angle_diff += Config::TWO_PI_F;
    }

    unwrapped_electrical_angle += angle_diff;
    previous_raw_electrical_angle = current_raw_electrical_angle;
}

void updateVelocityHistory(float velocity, unsigned long timestamp) {
    velocityHistory[velocityHistoryIndex] = {0.0f, velocity, timestamp};
    velocityHistoryIndex = (velocityHistoryIndex + 1) % Config::VELOCITY_HISTORY_SIZE;
    if (velocityHistoryIndex == 0) velocityHistoryFull = true;
}

void updateAngleHistory(float angle, unsigned long timestamp) {
    angleHistory[angleHistoryIndex] = {angle, 0.0f, timestamp};
    angleHistoryIndex = (angleHistoryIndex + 1) % Config::ANGLE_HISTORY_SIZE;
    if (angleHistoryIndex == 0) angleHistoryFull = true;
}

void updatePhaseHistory(float phaseA, float phaseB, float phaseC, unsigned long timestamp) {
    phaseHistory[phaseHistoryIndex] = {phaseA, phaseB, phaseC, timestamp};
    phaseHistoryIndex = (phaseHistoryIndex + 1) % Config::PHASE_HISTORY_SIZE;
    if (phaseHistoryIndex == 0) phaseHistoryFull = true;
}

float getSmoothedVelocity() {
    if (!velocityHistoryFull && velocityHistoryIndex < 3) return x_ekf(1); // Return EKF velocity if not enough history
    
    float sum_kahan = 0.0f;
    float compensation_kahan = 0.0f;
    float weight_sum = 0.0f;
    int count = velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex;
    
    // Apply exponential weighting (newer samples have higher weight) and Kahan summation
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float weight = exp(-i * 0.2f); // Exponential decay, tune 0.2 for decay rate
        kahanSum(velocityHistory[idx].velocity * weight, sum_kahan, compensation_kahan); // Use Kahan summation
        weight_sum += weight;
    }
    
    return (weight_sum > 0) ? sum_kahan / weight_sum : x_ekf(1); // Return EKF velocity if no valid sum
}

float getVelocityTrend() {
    if (!velocityHistoryFull && velocityHistoryIndex < 5) return 0.0f;
    
    int count = min(velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex, 5);
    float sum_x_kahan = 0.0f, comp_x = 0.0f;
    float sum_y_kahan = 0.0f, comp_y = 0.0f;
    float sum_xy_kahan = 0.0f, comp_xy = 0.0f;
    float sum_x2_kahan = 0.0f, comp_x2 = 0.0f;
    
    // Linear regression on recent velocity samples using Kahan summation
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float x = i; // Independent variable (time index)
        float y = velocityHistory[idx].velocity; // Dependent variable (velocity)
        kahanSum(x, sum_x_kahan, comp_x);
        kahanSum(y, sum_y_kahan, comp_y);
        kahanSum(x * y, sum_xy_kahan, comp_xy);
        kahanSum(x * x, sum_x2_kahan, comp_x2);
    }
    
    float denominator = count * sum_x2_kahan - sum_x_kahan * sum_x_kahan;
    return (fabsf(denominator) > 1e-6f) ? (count * sum_xy_kahan - sum_x_kahan * sum_y_kahan) / denominator : 0.0f;
}

bool isVelocityConsistent() {
    if (!velocityHistoryFull && velocityHistoryIndex < 3) return false;
    
    int count = min(velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex, 5);
    float variance_kahan = 0.0f;
    float comp_variance = 0.0f;
    float mean = smoothed_velocity; // Use the overall smoothed velocity as mean
    
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float diff = velocityHistory[idx].velocity - mean;
        kahanSum(diff * diff, variance_kahan, comp_variance); // Use Kahan summation for variance
    }
    variance_kahan /= count; // Calculate mean of squared differences
    
    return sqrt(variance_kahan) < Config::VELOCITY_NOISE_THRESHOLD;
}

// NEW FUNCTION: Infers rotation from phase history using linear regression
HistoricalPhaseStats getHistoricalPhaseStats() {
    HistoricalPhaseStats stats = {0.0f, 0.0f, 0.0f, false};

    if (!phaseHistoryFull && phaseHistoryIndex < 3) { // Need at least 3 points for a line
        return stats;
    }

    int count = phaseHistoryFull ? Config::PHASE_HISTORY_SIZE : phaseHistoryIndex;

    // Use relative timestamps (x values) and unwrapped angles (y values)
    float sum_x_kahan = 0.0f, comp_x = 0.0f;
    float sum_y_kahan = 0.0f, comp_y = 0.0f;
    float sum_xy_kahan = 0.0f, comp_xy = 0.0f;
    float sum_x2_kahan = 0.0f, comp_x2 = 0.0f;
    float sum_y2_kahan = 0.0f, comp_y2 = 0.0f;

    // To unwrap angles for regression, we need a reference point.
    float unwrapped_angles[Config::PHASE_HISTORY_SIZE];
    unsigned long relative_timestamps[Config::PHASE_HISTORY_SIZE];

    // Get the timestamp of the oldest entry in the current window to make timestamps relative
    unsigned long first_timestamp_in_window_idx = (phaseHistoryIndex - count + Config::PHASE_HISTORY_SIZE) % Config::PHASE_HISTORY_SIZE;
    unsigned long first_timestamp_in_window = phaseHistory[first_timestamp_in_window_idx].timestamp;
    
    // Initialize unwrapping with the first angle in the window
    float initial_raw_angle_unwrap;
    // Calculate the raw angle for the first entry in the window
    float a0 = phaseHistory[first_timestamp_in_window_idx].phaseA - Config::TARGET_QUIESCENT_AVG_ADC;
    float b0 = phaseHistory[first_timestamp_in_window_idx].phaseB - Config::TARGET_QUIESCENT_AVG_ADC;
    float c0 = phaseHistory[first_timestamp_in_window_idx].phaseC - Config::TARGET_QUIESCENT_AVG_ADC;
    float alpha0 = (2.0f / 3.0f) * (a0 - 0.5f * (b0 + c0));
    float beta0 = (2.0f / 3.0f) * ((b0 - c0) / Config::SQRT_3);
    initial_raw_angle_unwrap = atan2f(beta0, alpha0);
    unwrapped_angles[0] = initial_raw_angle_unwrap;
    relative_timestamps[0] = 0; // First point has relative time 0

    // Unroll angles and calculate sums for linear regression
    for (int i = 1; i < count; ++i) { // Start from second element
        int current_idx = (first_timestamp_in_window_idx + i) % Config::PHASE_HISTORY_SIZE;
        PhaseHistoryEntry current_entry = phaseHistory[current_idx];

        float a = current_entry.phaseA;
        float b = current_entry.phaseB;
        float c = current_entry.phaseC;
        
        a -= Config::TARGET_QUIESCENT_AVG_ADC;
        b -= Config::TARGET_QUIESCENT_AVG_ADC;
        c -= Config::TARGET_QUIESCENT_AVG_ADC;

        float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
        float beta = (2.0f / 3.0f) * ((b - c) / Config::SQRT_3);
        float raw_angle = atan2f(beta, alpha);

        // Unwrap angle relative to the previous unwrapped angle
        float angle_diff = raw_angle - unwrapped_angles[i-1];
        if (angle_diff > Config::PI_F) {
            angle_diff -= Config::TWO_PI_F;
        } else if (angle_diff < -Config::PI_F) {
            angle_diff += Config::TWO_PI_F;
        }
        unwrapped_angles[i] = unwrapped_angles[i-1] + angle_diff;
        
        relative_timestamps[i] = current_entry.timestamp - first_timestamp_in_window; // microsecond
    }

    // Now perform regression on (relative_timestamp_sec, unwrapped_angle)
    for (int i = 0; i < count; ++i) {
        float x = relative_timestamps[i] / 1e6f; // Convert to seconds
        float y = unwrapped_angles[i];

        kahanSum(x, sum_x_kahan, comp_x);
        kahanSum(y, sum_y_kahan, comp_y);
        kahanSum(x * y, sum_xy_kahan, comp_xy);
        kahanSum(x * x, sum_x2_kahan, comp_x2);
        kahanSum(y * y, sum_y2_kahan, comp_y2);
    }

    float N = static_cast<float>(count);
    float numerator = N * sum_xy_kahan - sum_x_kahan * sum_y_kahan;
    float denominator = N * sum_x2_kahan - sum_x_kahan * sum_x_kahan;

    if (fabsf(denominator) < 1e-9f) { // Avoid division by zero if all timestamps are the same or no variation
        return stats;
    }

    float slope = numerator / denominator; // This is the inferred angular velocity (rad/s)
    float intercept = (sum_y_kahan - slope * sum_x_kahan) / N; // Angle at t=0 of relative time scale

    // Calculate R-squared for consistency score
    float ss_tot_kahan = 0.0f, comp_ss_tot = 0.0f;
    float ss_res_kahan = 0.0f, comp_ss_res = 0.0f;
    float mean_y = sum_y_kahan / N;

    for (int i = 0; i < count; ++i) {
        float y_pred = slope * (relative_timestamps[i] / 1e6f) + intercept;
        kahanSum(powf(unwrapped_angles[i] - mean_y, 2), ss_tot_kahan, comp_ss_tot);
        kahanSum(powf(unwrapped_angles[i] - y_pred, 2), ss_res_kahan, comp_ss_res);
    }

    float r_squared = (ss_tot_kahan > 1e-9f) ? (1.0f - ss_res_kahan / ss_tot_kahan) : 1.0f; // Handle zero variance in y

    // Project angle to current time (micros())
    unsigned long current_relative_time_us = micros() - first_timestamp_in_window;
    float current_relative_time_sec = current_relative_time_us / 1e6f;
    float inferred_angle_at_current_time = slope * current_relative_time_sec + intercept;

    // Normalize this projected angle for consistent comparison, but regression was on unwrapped values
    stats.inferredAngleAtCurrentTime = normalizeAngle(inferred_angle_at_current_time);
    stats.inferredVelocity = slope;
    stats.consistencyScore = constrain(r_squared, 0.0f, 1.0f); // R-squared can sometimes be slightly negative if model is worse than mean
    stats.isValid = true;

    return stats;
}


bool isRotationStopping() {
    // Check if velocity is decreasing consistently
    float trend = getVelocityTrend();
    float current_vel = fabsf(smoothed_velocity);
    
    // Rotation is stopping if:
    // 1. Velocity magnitude is small and decreasing
    // 2. Velocity trend is negative (slowing down)
    // 3. Velocity is consistent (not noisy)
    
    return (current_vel < Config::MIN_VELOCITY_CONSISTENCY && 
            trend < -0.01f && 
            isVelocityConsistent());
}

void applyEKF(float measured_angle) {
    unsigned long now = micros();
    float dt = (now - prev_ekf_micros) / 1e6f; // Time step in seconds
    // Ensure dt is not zero to prevent division by zero in some cases
    if (dt == 0) { 
        dt = 1e-6f; // Assign a small non-zero value if micros() hasn't updated
    }
    prev_ekf_micros = now;

    float angle = x_ekf(0), velocity = x_ekf(1);

    // Predict
    float predicted_angle = normalizeAngle(angle + velocity * dt);
    float predicted_velocity = velocity; // Assuming constant velocity for prediction

    BLA::Matrix<N_STATES_EKF, 1> x_pred;
    x_pred(0) = predicted_angle;
    x_pred(1) = predicted_velocity;

    // State Transition Matrix (A)
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A;
    A(0, 0) = 1; A(0, 1) = dt;
    A(1, 0) = 0; A(1, 1) = 1;

    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A_transposed = A.Transpose();

    P_ekf = A * P_ekf * A_transposed + Q_ekf; // Predicted covariance

    // --- Novel Method Integration: Use Historical Phase Stats for R_ekf Adjustment ---
    HistoricalPhaseStats historicalStats = getHistoricalPhaseStats();
    historical_inferred_velocity = historicalStats.inferredVelocity;
    historical_inferred_angle = historicalStats.inferredAngleAtCurrentTime;
    phase_history_consistency_score = historicalStats.consistencyScore;
    
    // Determine phase consistency based on a threshold of R-squared
    is_consistent_phases = historicalStats.isValid && (historical_history_consistency_score > Config::PHASE_CONSISTENCY_R_SQUARED_THRESHOLD);

    // Dynamically adjust R_ekf based on phase history consistency.
    // When consistency is high (score near 1), R_ekf is small (trust measurement).
    // When consistency is low (score near 0), R_ekf is large (trust prediction more).
    // Interpolate between R_ANGLE_EKF (min R) and R_ANGLE_EKF_MAX (max R).
    R_ekf(0, 0) = Config::R_ANGLE_EKF + (1.0f - phase_history_consistency_score) * (Config::R_ANGLE_EKF_MAX - Config::R_ANGLE_EKF);
    
    // If there's not enough valid history for regression, set R_ekf to its max to favor prediction
    if (!historicalStats.isValid) { 
        R_ekf(0, 0) = Config::R_ANGLE_EKF_MAX;
    }

    // Update Step
    float y = normalizeAngle(measured_angle - predicted_angle); // Measurement residual (innovation)
    BLA::Matrix<1, 1> S = R_ekf + P_ekf(0, 0); // Innovation covariance (H is [1 0], so H*P*H' is just P(0,0))
    BLA::Matrix<N_STATES_EKF, 1> K; // Kalman Gain
    K(0) = P_ekf(0, 0) / S(0, 0);
    K(1) = P_ekf(1, 0) / S(0, 0);
    
    x_ekf = x_pred + K * y; // Updated state estimate
    x_ekf(0) = normalizeAngle(x_ekf(0)); // Keep angle normalized

    // Update covariance matrix
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_KH = I_matrix_ekf; // (I - KH)
    I_KH(0, 0) -= K(0, 0); // H is implicitly [1 0] for angle measurement
    I_KH(1, 0) -= K(1, 0);
    P_ekf = I_KH * P_ekf; // Updated covariance

    // Update history for smoothed velocity and angle calculation
    updateVelocityHistory(x_ekf(1), now);
    updateAngleHistory(x_ekf(0), now); // EKF's wrapped angle history
    
    // Update diagnostic flags
    smoothed_velocity = getSmoothedVelocity();
    velocity_trend = getVelocityTrend();
    is_stopping = isRotationStopping(); 
    bool is_consistent_vel = isVelocityConsistent(); 
    float vel_magnitude = fabsf(x_ekf(1)); // Use EKF's current velocity for magnitude check
    bool is_moving = vel_magnitude > Config::VELOCITY_STEP_THRESHOLD;

    // Update EKF lock status (for confidence indication only)
    // The lock criteria now heavily depends on the phase_history_consistency_score
    if (is_moving && is_consistent_vel && is_consistent_phases && !is_stopping) {
        if (!ekfLocked) { // If not yet locked
            consistent_direction_count++; // Increment counter for consistent movement
            if (consistent_direction_count > 3) { // Require a few consistent readings
                ekfLocked = true;
                Serial.println(F("EKF locked - confident movement detected."));
            }
        }
        lastNonZeroVelMicros = now;
    } else { // Conditions for unlocking EKF
        if (ekfLocked) {
            bool should_unlock = false;
            if (is_stopping) {
                should_unlock = true;
                Serial.println(F("EKF unlocked - rotation stopping detected."));
            } else if (!is_moving && (now - lastNonZeroVelMicros > Config::ZERO_VEL_LOCK_TIMEOUT)) {
                should_unlock = true;
                Serial.println(F("EKF unlocked - prolonged zero velocity."));
            } else if (!is_consistent_vel && vel_magnitude < Config::MIN_VELOCITY_CONSISTENCY) {
                should_unlock = true;
                Serial.println(F("EKF unlocked - inconsistent low velocity."));
            } else if (!is_consistent_phases) { // Now relies on the R-squared threshold
                should_unlock = true;
                Serial.println(F("EKF unlocked - inconsistent phase signals (low R-squared)."));
            }
            if (should_unlock) {
                ekfLocked = false;
                consistent_direction_count = 0; // Reset counter upon unlock
            }
        } else {
             consistent_direction_count = 0; // Reset if not locked and not meeting criteria
        }
    }
}

// Normalizes an angle to the range (-PI, PI]
float normalizeAngle(float a) {
    while (a <= -Config::PI_F) a += Config::TWO_PI_F;
    while (a > Config::PI_F) a -= Config::TWO_PI_F;
    return a;
}

// Safe analogRead to prevent negative or overflow values
int checkedAnalogRead(int pin) {
    int val = analogRead(pin);
    if (val < 0) val = 0;
    if (val > 1023) val = 1023;
    return val;
}

void logData() {
    // Serial Plotter Output: 
    // 0: KF_PhaseA_Adjusted, 
    // 1: KF_PhaseB_Adjusted, 
    // 2: KF_PhaseC_Adjusted, 
    // 3: Unwrapped_Angle (EKF's angle, unwrapped), 
    // 4: EKF_Velocity (EKF's direct velocity), 
    // 5: EKF_Locked_Status, 
    // 6: Phase_Consistency_Status (from R-squared), 
    // 7: Estimated_Mechanical_Steps, 
    // 8: Revolution_Count,
    // 9: Inferred_Angle_History (Candidate 2 Angle from Phase History Regression),
    // 10: Inferred_Velocity_History (Candidate 2 Velocity from Phase History Regression)
    Serial.print(x_hat_kf1(0)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(x_hat_kf1(1)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(x_hat_kf1(2)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(unwrapped_electrical_angle, 3); Serial.print(","); 
    Serial.print(x_ekf(1), 4); Serial.print(","); 
    Serial.print(ekfLocked ? 1 : 0); Serial.print(","); 
    Serial.print(is_consistent_phases ? 1 : 0);Serial.print(","); 
    Serial.print(estimated_mechanical_steps);Serial.print(",");
    Serial.print(revolution_count);Serial.print(",");
    Serial.print(historical_inferred_angle, 3); Serial.print(","); 
    Serial.println(historical_inferred_velocity, 4); 

    // Uncomment the lines below if you also want detailed text output in Serial Monitor:
    /*
    Serial.print("KF_A: "); Serial.print(x_hat_kf1(0), 2);
    Serial.print("\tKF_B: "); Serial.print(x_hat_kf1(1), 2);
    Serial.print("\tKF_C: "); Serial.print(x_hat_kf1(2), 2);
    Serial.print("\tUnwrapped Ang: "); Serial.print(unwrapped_electrical_angle, 3);
    Serial.print("\tEKF Vel: "); Serial.print(x_ekf(1), 4);
    Serial.print("\tHist Ang: "); Serial.print(historical_inferred_angle, 3);
    Serial.print("\tHist Vel: "); Serial.print(historical_inferred_velocity, 4);
    Serial.print("\tPhase Consist: "); Serial.print(phase_history_consistency_score, 2);
    Serial.print("\tEKF Lock: "); Serial.print(ekfLocked ? "Locked" : "Unlocked");
    Serial.print("\tIs Consistent Phases: "); Serial.print(is_consistent_phases ? "Yes" : "No");
    Serial.print("\tMech Steps: "); Serial.print(estimated_mechanical_steps);
    Serial.print("\tRevs: "); Serial.println(revolution_count);
    */
}
