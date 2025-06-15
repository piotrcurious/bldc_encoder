#include <Arduino.h>
#include <BasicLinearAlgebra.h>

namespace Config {
    constexpr int PHASE_A_PIN = A0, PHASE_B_PIN = A1, PHASE_C_PIN = A2, COMMON_PIN = 9;
    constexpr float ELECTRICAL_STEPS_PER_REV = 1.0f;
    // Define PI and TWO_PI_F for clarity and consistency
    const float PI_F = 3.14159265358979323846f;
    constexpr float TWO_PI_F = 2.0f * PI; // PI is defined in arduino framework
    constexpr float SQRT_3 = 1.7320508f;
    
    // Original NUM_SAMPLES is now a default/initial value
    constexpr int DEFAULT_NUM_SAMPLES = 16; 

    // New parameters for dynamic oversampling
    constexpr int MIN_NUM_SAMPLES = 16;   // Minimum number of samples
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
    constexpr float Q_ANGLE_EKF = 1.1f, Q_VEL_EKF = 80.1f, R_ANGLE_EKF = 1.1f; // Adjusted noise parameters
    
    // New parameters for improved tracking
    constexpr int VELOCITY_HISTORY_SIZE = 4;
    constexpr int ANGLE_HISTORY_SIZE = 4;
    constexpr int PHASE_HISTORY_SIZE = 4; // History for KF filtered phases
    constexpr float MIN_VELOCITY_CONSISTENCY = 0.0003f; // Minimum velocity for consistent direction
    constexpr float VELOCITY_NOISE_THRESHOLD = 0.006f; // Threshold for velocity noise rejection
    constexpr float ANGLE_HISTORY_NOISE_THRESHOLD = 0.0005f; // Threshold for angle consistency from phase history (radians)
    constexpr unsigned long MIN_STABLE_TIME_US = 500UL; // Minimum time for stable velocity
}

#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1;

#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf;
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf;
BLA::Matrix<1, 1> R_ekf;

// Enhanced state tracking
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

HistoryEntry velocityHistory[Config::VELOCITY_HISTORY_SIZE];
HistoryEntry angleHistory[Config::ANGLE_HISTORY_SIZE];
PhaseHistoryEntry phaseHistory[Config::PHASE_HISTORY_SIZE]; // New phase history buffer
int velocityHistoryIndex = 0, angleHistoryIndex = 0, phaseHistoryIndex = 0; // New index for phase history
bool velocityHistoryFull = false, angleHistoryFull = false, phaseHistoryFull = false; // New flag for phase history

// Global variables for the new rotation tracking system
float previous_raw_electrical_angle = 0.0f; // Stores the previous raw angle from KF for unwrapping
bool first_angle_reading = true; // Flag to initialize unwrapped angle on first read

// Main angle tracking variables
float unwrapped_electrical_angle = 0.0f; // This is the continuous, unbounded angle
long revolution_count = 0;
long estimated_mechanical_steps = 0;

bool ekfLocked = false; // EKF lock indicates confidence
unsigned long lastNonZeroVelMicros = 0, prev_ekf_micros = 0;

float smoothed_velocity = 0.0f;
float velocity_trend = 0.0f;
unsigned long last_direction_change = 0;
int consistent_direction_count = 0;
float last_stable_velocity = 0.0f;
bool is_stopping = false;
bool is_consistent_phases = false; // Flag to indicate phase consistency

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
bool checkPhaseConsistency();
bool isVelocityConsistent();
bool isRotationStopping();

void updateUnwrappedElectricalAngle(float current_raw_electrical_angle);

int getDynamicNumSamples(float current_velocity);

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

    // 3. Update the continuous, unwrapped electrical angle
    updateUnwrappedElectricalAngle(current_raw_angle_from_kf);

    // 4. Apply EKF: EKF now filters the current wrapped angle and estimates velocity
    applyEKF(current_raw_angle_from_kf); 

    // 5. Calculate mechanical steps from the unwrapped angle
    // Each electrical revolution (2*PI) corresponds to ELECTRICAL_STEPS_PER_REV mechanical steps.
    // E.g., for a 3-phase, 2-pole pair motor, 1 electrical rev = 1 mechanical rev.
    // For 3-phase, 3-pole pair motor, 1 electrical rev = 2/3 mechanical rev.
    // So mechanical revolutions = unwrapped_electrical_angle / Config::TWO_PI_F
    // And mechanical steps = mechanical revolutions * Config::ELECTRICAL_STEPS_PER_REV
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
                pinMode(Config::COMMON_PIN, INPUT);
                delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
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
            is_consistent_phases = false;
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
            pinMode(Config::COMMON_PIN, INPUT);
            delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
            readings_tune[0] += checkedAnalogRead(Config::PHASE_A_PIN);
            readings_tune[1] += checkedAnalogRead(Config::PHASE_B_PIN);
            readings_tune[2] += checkedAnalogRead(Config::PHASE_C_PIN);
        }
        for (int i = 0; i < 3; i++) readings_tune[i] /= auto_tune_samples;
        applyMultivariateKalmanFilter(readings_tune);
        
        int current_avg_adc = static_cast<int>((x_hat_kf1(0) + x_hat_kf1(1) + x_hat_kf1(2)) / 3.0f);
        int error = current_avg_adc - Config::TARGET_QUIESCENT_AVG_ADC;

       // if (iteration % (Config::PWM_CALIBRATION_ITERATIONS / 10 + 1) == 0) { // Reduce serial output frequency
        //    Serial.print(F("Step ")); Serial.print(iteration);
        //    Serial.print(F(" | PWM: ")); Serial.print(Config::EXCITATION_PWM_VALUE);
        //    Serial.print(F(" | Avg ADC: ")); Serial.print(current_avg_adc);
        //    Serial.print(F(" | Error: ")); Serial.println(error);
        //}

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

// **MODIFIED** Function for dynamic oversampling - now takes velocity
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
    
    // **MODIFIED:** Use the EKF's estimated velocity for dynamic oversampling
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

        // Step 5: Stop PWM and make COMMON_PIN float (high impedance).
        analogWrite(Config::COMMON_PIN, 0); // Stop PWM output
        //pinMode(Config::COMMON_PIN, INPUT); // Make the pin float (high impedance)

        // Step 6: Wait for settling time after the common pin is floating.
        //delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);

        // Step 7: Perform ADC readings while phases are floating.
        readings[0] += checkedAnalogRead(Config::PHASE_A_PIN);
        readings[1] += checkedAnalogRead(Config::PHASE_B_PIN);
        readings[2] += checkedAnalogRead(Config::PHASE_C_PIN);
    }
    // Average the readings over the number of samples
    for (int i = 0; i < 3; i++) readings[i] /= num_samples_to_take;
    // Apply the Kalman Filter to the averaged readings
    applyMultivariateKalmanFilter(readings);
        pinMode(Config::COMMON_PIN, OUTPUT);
        digitalWrite(Config::COMMON_PIN, LOW); // set up common pin low for rest of the time 

}

void applyMultivariateKalmanFilter(const long readings[3]) {
    BLA::Matrix<N_STATES_KF1, 1> z, y;
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S; // S must be N_STATES_KF1 x N_STATES_KF1 (3x3)

    for (int i = 0; i < 3; i++) z(i) = readings[i];
    y = z - x_hat_kf1;
    S = P_kf1 + R_kf1; // Now S is correctly 3x3
    auto K = P_kf1 * BLA::Inverse(S); // Inverse will now work on 3x3 S
    x_hat_kf1 += K * y;
    P_kf1 = (I_matrix_kf1 - K) * P_kf1 + Q_kf1;

    // Update KF phase history
    updatePhaseHistory(x_hat_kf1(0), x_hat_kf1(1), x_hat_kf1(2), micros());
}

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

// **NEW** Function for continuous angle unwrapping
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
    if (!velocityHistoryFull && velocityHistoryIndex < 3) return x_ekf(1);
    
    float sum_kahan = 0.0f;
    float compensation_kahan = 0.0f;
    float weight_sum = 0.0f;
    int count = velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex;
    
    // Apply exponential weighting (newer samples have higher weight) and Kahan summation
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float weight = exp(-i * 0.2f); // Exponential decay
        kahanSum(velocityHistory[idx].velocity * weight, sum_kahan, compensation_kahan); // Use Kahan summation
        weight_sum += weight;
    }
    
    return (weight_sum > 0) ? sum_kahan / weight_sum : x_ekf(1);
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
        float x = i;
        float y = velocityHistory[idx].velocity;
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
    float mean = smoothed_velocity;
    
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float diff = velocityHistory[idx].velocity - mean;
        kahanSum(diff * diff, variance_kahan, comp_variance); // Use Kahan summation for variance
    }
    variance_kahan /= count;
    
    return sqrt(variance_kahan) < Config::VELOCITY_NOISE_THRESHOLD;
}

bool checkPhaseConsistency() {
    if (!phaseHistoryFull && phaseHistoryIndex < 5) return false; // Need some history to check consistency
    
    int count = phaseHistoryFull ? Config::PHASE_HISTORY_SIZE : phaseHistoryIndex;
    
    // Calculate angles from history and check standard deviation
    float sum_angles_kahan = 0.0f;
    float comp_angles = 0.0f;
    float angles[Config::PHASE_HISTORY_SIZE]; // Use max size for array, actual elements used by 'count'

    for (int i = 0; i < count; ++i) {
        int idx = (phaseHistoryIndex - 1 - i + Config::PHASE_HISTORY_SIZE) % Config::PHASE_HISTORY_SIZE;
        float a = phaseHistory[idx].phaseA;
        float b = phaseHistory[idx].phaseB;
        float c = phaseHistory[idx].phaseC;
        
        // Calculate inferred angle from historical KF outputs
        // IMPORTANT: Subtract the ADC offset here as well for consistency with getInferredElectricalAngle
        a -= Config::TARGET_QUIESCENT_AVG_ADC;
        b -= Config::TARGET_QUIESCENT_AVG_ADC;
        c -= Config::TARGET_QUIESCENT_AVG_ADC;

        float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
        float beta = (2.0f / 3.0f) * ((b - c) / Config::SQRT_3);
        angles[i] = atan2f(beta, alpha);
        kahanSum(angles[i], sum_angles_kahan, comp_angles); // Use Kahan summation
    }

    float mean_angle = sum_angles_kahan / count; // Use Kahan sum for mean
    float variance_kahan = 0.0f;
    float comp_variance = 0.0f;
    for (int i = 0; i < count; ++i) {
        float diff = normalizeAngle(angles[i] - mean_angle); // Normalize difference for circular data
        kahanSum(diff * diff, variance_kahan, comp_variance); // Use Kahan summation for variance
    }
    variance_kahan /= count;
    
    return sqrt(variance_kahan) < Config::ANGLE_HISTORY_NOISE_THRESHOLD;
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
    prev_ekf_micros = now;

    float angle = x_ekf(0), velocity = x_ekf(1);

    // Predict
    float predicted_angle = normalizeAngle(angle + velocity * dt);
    float predicted_velocity = velocity;

    BLA::Matrix<N_STATES_EKF, 1> x_pred;
    x_pred(0) = predicted_angle;
    x_pred(1) = predicted_velocity;

    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A;
    A(0, 0) = 1; A(0, 1) = dt;
    A(1, 0) = 0; A(1, 1) = 1;

    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A_transposed;
    A_transposed(0, 0) = A(0, 0); A_transposed(0, 1) = A(1, 0);
    A_transposed(1, 0) = A(0, 1); A_transposed(1, 1) = A(1, 1);

    P_ekf = A * P_ekf * A_transposed + Q_ekf;

    // Update
    float y = normalizeAngle(measured_angle - predicted_angle); // Measurement residual
    BLA::Matrix<1, 1> S = R_ekf + P_ekf(0, 0); // Innovation covariance
    BLA::Matrix<N_STATES_EKF, 1> K; // Kalman Gain
    K(0) = P_ekf(0, 0) / S(0, 0);
    K(1) = P_ekf(1, 0) / S(0, 0);
    
    x_ekf = x_pred + K * y; // Updated state estimate

    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_KH; // (I - KH)
    I_KH = I_matrix_ekf;
    I_KH(0, 0) -= K(0, 0);
    I_KH(1, 0) -= K(1, 0);
    P_ekf = I_KH * P_ekf; // Updated covariance

    // Update history for smoothed velocity calculation
    updateVelocityHistory(x_ekf(1), now);
    updateAngleHistory(x_ekf(0), now); // EKF's wrapped angle history
    
    // Update diagnostic flags (these no longer control unwrapped_electrical_angle directly)
    smoothed_velocity = getSmoothedVelocity();
    velocity_trend = getVelocityTrend();
    is_stopping = isRotationStopping();
    is_consistent_phases = checkPhaseConsistency();
    bool is_consistent_vel = isVelocityConsistent();
    float vel_magnitude = fabsf(smoothed_velocity);
    bool is_moving = vel_magnitude > Config::VELOCITY_STEP_THRESHOLD;

    // Update EKF lock status (for confidence indication only)
    if (is_moving && is_consistent_vel && is_consistent_phases && !is_stopping) {
        if (!ekfLocked && consistent_direction_count > 3) {
            ekfLocked = true;
            Serial.println(F("EKF locked - confident movement detected."));
        }
        lastNonZeroVelMicros = now;
    } else if (ekfLocked) {
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
        } else if (!is_consistent_phases) {
            should_unlock = true;
            Serial.println(F("EKF unlocked - inconsistent phase signals."));
        }
        if (should_unlock) {
            ekfLocked = false;
            consistent_direction_count = 0;
        }
    }
}

// Normalizes an angle to the range (-PI, PI]
float normalizeAngle(float a) {
    while (a <= -Config::PI_F) a += Config::TWO_PI_F;
    while (a > Config::PI_F) a -= Config::TWO_PI_F;
    return a;
}

int checkedAnalogRead(int pin) {
    int val = analogRead(pin);
    if (val < 0) val = 0;
    if (val > 1023) val = 1023;
    return val;
}

void logData() {
    // Serial Plotter Output: KF_PhaseA, KF_PhaseB, KF_PhaseC, Unwrapped_Angle, Smoothed_Velocity, EKF_Locked_Status, Phase_Consistency_Status, Estimated_Mechanical_Steps, Revolution_Count
    Serial.print(x_hat_kf1(0)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(x_hat_kf1(1)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(x_hat_kf1(2)-Config::TARGET_QUIESCENT_AVG_ADC); Serial.print(",");
    Serial.print(unwrapped_electrical_angle, 3); Serial.print(","); // Primary continuous angle
    Serial.print(smoothed_velocity, 4); Serial.print(","); // Smoothed Velocity from EKF
    Serial.print(ekfLocked ? 1 : 0); Serial.print(","); // EKF Lock status
    Serial.print(is_consistent_phases ? 1 : 0);Serial.print(","); // Phase Consistency status
    Serial.print(estimated_mechanical_steps);Serial.print(",");
    Serial.println(revolution_count);

    // Uncomment the lines below if you also want detailed text output in Serial Monitor:
    /*
    Serial.print("KF_A: "); Serial.print(x_hat_kf1(0), 2);
    Serial.print("\tKF_B: "); Serial.print(x_hat_kf1(1), 2);
    Serial.print("\tKF_C: "); Serial.print(x_hat_kf1(2), 2);
    Serial.print("\tUnwrapped Ang: "); Serial.print(unwrapped_electrical_angle, 3);
    Serial.print("\tSmoothed Vel: "); Serial.print(smoothed_velocity, 4);
    Serial.print("\tEKF Lock: "); Serial.print(ekfLocked ? "Locked" : "Unlocked");
    Serial.print("\tPhase Consist: "); Serial.print(is_consistent_phases ? "Consistent" : "Inconsistent");
    Serial.print("\tMech Steps: "); Serial.print(estimated_mechanical_steps);
    Serial.print("\tRevs: "); Serial.println(revolution_count);
    */
}
