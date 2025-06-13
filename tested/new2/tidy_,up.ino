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
    constexpr int NUM_SAMPLES_PER_PHASE = 16; // Number of samples for averaging each phase
    int EXCITATION_PWM_VALUE = 1; // Calibrated during auto-tuning, mutable
    constexpr int EXCITATION_PULSE_WIDTH_US = 1000; // Duration for excitation pulse
    constexpr int FLOATING_SETTLING_TIME_US = 10;   // Time to settle after floating

    // Autotuning Parameters
    constexpr int TARGET_QUIESCENT_AVG_ADC = 80;
    constexpr int PWM_ADJUST_STEP = 1;
    constexpr int PWM_CALIBRATION_ITERATIONS = 2000;
    constexpr int PWM_CALIBRATION_TOLERANCE = 1;

    // EKF & KF Noise Parameters (Adjusted for potentially better performance)
    constexpr float Q_DIAGONAL_VALUE_KF1 = 0.5f;   // Process noise for KF1 (lower for less trust in process)
    constexpr float R_DIAGONAL_VALUE_KF1 = 200.0f; // Measurement noise for KF1 (higher for more noise in ADC)

    constexpr float Q_ANGLE_EKF = 0.0005f; // Process noise for angle (lower for smoother angle)
    constexpr float Q_VEL_EKF = 0.005f;    // Process noise for velocity (lower for smoother velocity)
    constexpr float R_ANGLE_EKF = 0.2f;    // Measurement noise for EKF (higher for more noise in inferred angle)

    // Velocity & Locking Logic Parameters
    constexpr int VELOCITY_HISTORY_SIZE = 15; // Increased history for better smoothing/trend detection
    constexpr int ANGLE_HISTORY_SIZE = 7;     // Increased history for better consistency checks
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
}

// --- Global Variables ---
#define N_STATES_KF1 3
BLA::Matrix<N_STATES_KF1, 1> x_hat_kf1;                  // KF1 state: [Phase A, Phase B, Phase C]
BLA::Matrix<N_STATES_KF1, N_STATES_KF1> P_kf1, Q_kf1, R_kf1, I_matrix_kf1;

#define N_STATES_EKF 2
BLA::Matrix<N_STATES_EKF, 1> x_ekf;                      // EKF state: [Angle, Velocity]
BLA::Matrix<N_STATES_EKF, N_STATES_EKF> P_ekf, Q_ekf, I_matrix_ekf;
BLA::Matrix<1, 1> R_ekf;

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
void executeSensingCycle();
void applyMultivariateKalmanFilter(const int readings[3]);
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
    // For ATmega328P: ADPS2=1, ADPS1=0, ADPS0=0 --> Prescaler = 16 (1MHz ADC clock)
    // For ATmega328P: ADPS2=1, ADPS1=0, ADPS0=1 --> Prescaler = 32 (500kHz ADC clock)
    ADCSRA = (ADCSRA & ~0x07) | 0x05; // Set prescaler to 32 (0x05 for 32, 0x04 for 16)
    // Configure Timer1 for high-speed PWM on pin 9 (COMMON_PIN)
    // This assumes COMMON_PIN is associated with Timer1, e.g., OC1A on Uno.
    // For Fast PWM, 8-bit mode (WGM12, WGM10)
    // TCCR1A = _BV(COM1A1) | _BV(WGM10); // Clear OC1A on Compare Match, set on BOTTOM, Fast PWM 8-bit
    // TCCR1B = _BV(WGM12) | _BV(CS10);    // Fast PWM 8-bit, No prescaling (full speed)
    // The existing TCCR1B = (TCCR1B & 0b11111000) | 0x01; sets Timer1 prescaler to 1, which is fine for direct PWM.
    // Ensure pinMode(COMMON_PIN, OUTPUT); is done.
    pinMode(Config::COMMON_PIN, OUTPUT);
    digitalWrite(Config::COMMON_PIN, LOW); // Start with common pin low

    autoTuneExcitationPwm(); // Calibrate excitation PWM

    // Initialize KF1 state with the calibrated quiescent ADC value
    x_hat_kf1.Fill((float)Config::TARGET_QUIESCENT_AVG_ADC);
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
    x_ekf(1) = 0.0f;                         // Initial velocity is zero
    P_ekf.Fill(0.0f);
    P_ekf(0, 0) = 1.0f; // Initial angle uncertainty
    P_ekf(1, 1) = 1.0f; // Initial velocity uncertainty
    prev_ekf_micros = micros();
    unwrapped_electrical_angle = x_ekf(0);

    Serial.println(F("Enhanced system ready. Rotate motor to begin estimation."));
}

// --- Main Loop Function ---
void loop() {
    executeSensingCycle();
    float inferred_angle = getInferredElectricalAngle();
    applyEKF(inferred_angle);

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
            x_ekf(1) = 0.0f;                         // Reset EKF velocity
            P_ekf(0, 0) = 1.0f; P_ekf(1, 1) = 1.0f;  // Reset EKF uncertainty
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

void autoTuneExcitationPwm() {
    int sum = 0;
    // No need for prev_avg, direct comparison with TARGET_QUIESCENT_AVG_ADC is sufficient
    Serial.print(F("Auto-tuning PWM. Target ADC: ")); Serial.println(Config::TARGET_QUIESCENT_AVG_ADC);

    for (int iter = 0; iter < Config::PWM_CALIBRATION_ITERATIONS; iter++) {
        sum = 0;
        // Collect multiple samples for better average
        for (int i = 0; i < Config::NUM_SAMPLES_PER_PHASE; i++) {
            // Set common pin to output LOW before analogWrite for consistency
            digitalWrite(Config::COMMON_PIN, LOW); // Ensures common pin is initially low
            analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
            delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);
            // Switch common pin to input (floating) for measurement
            pinMode(Config::COMMON_PIN, INPUT);
            delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);
            int reading = checkedAnalogRead(Config::PHASE_A_PIN); // Measure one phase
            sum += reading;
        }
        int avg = sum / Config::NUM_SAMPLES_PER_PHASE;

        if (abs(avg - Config::TARGET_QUIESCENT_AVG_ADC) <= Config::PWM_CALIBRATION_TOLERANCE) {
            Serial.print(F("PWM calibration converged in ")); Serial.print(iter + 1); Serial.println(F(" iterations."));
            break; // Converged
        }

        // Adjust PWM value
        if (avg < Config::TARGET_QUIESCENT_AVG_ADC) {
            Config::EXCITATION_PWM_VALUE += Config::PWM_ADJUST_STEP;
        } else {
            Config::EXCITATION_PWM_VALUE -= Config::PWM_ADJUST_STEP;
        }
        Config::EXCITATION_PWM_VALUE = constrain(Config::EXCITATION_PWM_VALUE, 0, 255);

        // Optional: print progress during tuning
        // if (iter % 100 == 0) {
        //     Serial.print(F("Iter: ")); Serial.print(iter);
        //     Serial.print(F(", Avg ADC: ")); Serial.print(avg);
        //     Serial.print(F(", PWM: ")); Serial.println(Config::EXCITATION_PWM_VALUE);
        // }
    }
    Serial.print(F("Calibrated PWM value: "));
    Serial.println(Config::EXCITATION_PWM_VALUE);
    // Ensure COMMON_PIN is returned to INPUT state after calibration
    pinMode(Config::COMMON_PIN, INPUT);
}

void executeSensingCycle() {
    int readings[3] = {0}; // Raw ADC sums for each phase
    int phasePins[] = {Config::PHASE_A_PIN, Config::PHASE_B_PIN, Config::PHASE_C_PIN};

    for (int i = 0; i < Config::NUM_SAMPLES_PER_PHASE; i++) {
        digitalWrite(Config::COMMON_PIN, LOW); // Drive common low for excitation
        analogWrite(Config::COMMON_PIN, Config::EXCITATION_PWM_VALUE);
        delayMicroseconds(Config::EXCITATION_PULSE_WIDTH_US);

        pinMode(Config::COMMON_PIN, INPUT); // Set common to floating
        delayMicroseconds(Config::FLOATING_SETTLING_TIME_US);

        // Read all three phases while common is floating
        readings[0] += checkedAnalogRead(phasePins[0]);
        readings[1] += checkedAnalogRead(phasePins[1]);
        readings[2] += checkedAnalogRead(phasePins[2]);
    }

    // Average the readings
    for (int i = 0; i < 3; i++) {
        readings[i] /= Config::NUM_SAMPLES_PER_PHASE;
    }

    applyMultivariateKalmanFilter(readings);
}

void applyMultivariateKalmanFilter(const int readings[3]) {
    BLA::Matrix<N_STATES_KF1, 1> z; // Measurement vector
    for (int i = 0; i < 3; i++) z(i) = readings[i];

    // Kalman Filter equations
    BLA::Matrix<N_STATES_KF1, 1> y = z - x_hat_kf1; // Innovation
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> S = P_kf1 + R_kf1; // Innovation covariance
    BLA::Matrix<N_STATES_KF1, N_STATES_KF1> K = P_kf1 * BLA::Inverse(S); // Kalman gain

    x_hat_kf1 += K * y; // Update state estimate
    P_kf1 = (I_matrix_kf1 - K) * P_kf1 + Q_kf1; // Update error covariance
}

float getInferredElectricalAngle() {
    // Clarke transform for alpha-beta components
    float a = x_hat_kf1(0);
    float b = x_hat_kf1(1);
    float c = x_hat_kf1(2);

    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta = (1.0f / Config::SQRT_3) * (b - c); // Simplified calculation (2/3 * sqrt(3) / 2)

    return atan2f(beta, alpha); // Returns angle in (-PI, PI]
}

void updateVelocityHistory(float velocity, unsigned long timestamp) {
    velocityHistory[velocityHistoryIndex] = {0.0f, velocity, timestamp}; // Angle field unused here
    velocityHistoryIndex = (velocityHistoryIndex + 1) % Config::VELOCITY_HISTORY_SIZE;
    if (velocityHistoryIndex == 0) velocityHistoryFull = true;
}

void updateAngleHistory(float angle, unsigned long timestamp) {
    angleHistory[angleHistoryIndex] = {angle, 0.0f, timestamp}; // Velocity field unused here
    angleHistoryIndex = (angleHistoryIndex + 1) % Config::ANGLE_HISTORY_SIZE;
    if (angleHistoryIndex == 0) angleHistoryFull = true;
}

float getSmoothedVelocity() {
    if (!velocityHistoryFull && velocityHistoryIndex < 3) { // Need at least 3 samples to smooth
        return x_ekf(1); // Return current EKF velocity if not enough history
    }

    float sum = 0.0f;
    float weight_sum = 0.0f;
    int count = velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex;

    // Apply exponential weighting (newer samples have higher weight)
    for (int i = 0; i < count; i++) {
        // Calculate index for reverse chronological order
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float weight = exp(-i * 0.2f); // Exponential decay factor, adjust as needed
        sum += velocityHistory[idx].velocity * weight;
        weight_sum += weight;
    }

    return (weight_sum > 0) ? sum / weight_sum : x_ekf(1);
}

float getVelocityTrend() {
    // Need a minimum of 5 samples for a meaningful linear regression
    if (!velocityHistoryFull && velocityHistoryIndex < 5) {
        return 0.0f;
    }

    int count = min(velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex, 5); // Use last 5 samples
    float sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;

    // Linear regression on recent velocity samples
    // x represents the time index (0 for newest, 1 for second newest, etc.)
    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float x = (float)i; // Independent variable (time/sample index)
        float y = velocityHistory[idx].velocity; // Dependent variable (velocity)
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
    }

    float denominator = count * sum_x2 - sum_x * sum_x;
    // Calculate slope (trend)
    return (fabsf(denominator) > 1e-6f) ? (count * sum_xy - sum_x * sum_y) / denominator : 0.0f;
}

// Check if the current velocity is consistent (low noise)
bool isVelocityConsistent(float velocity_to_check) {
    if (!velocityHistoryFull && velocityHistoryIndex < 3) return false;

    int count = min(velocityHistoryFull ? Config::VELOCITY_HISTORY_SIZE : velocityHistoryIndex, 5); // Check last 5 samples
    float variance = 0.0f;
    // Use the `velocity_to_check` as the mean for variance calculation, often `smoothed_velocity`
    float mean = velocity_to_check;

    for (int i = 0; i < count; i++) {
        int idx = (velocityHistoryIndex - 1 - i + Config::VELOCITY_HISTORY_SIZE) % Config::VELOCITY_HISTORY_SIZE;
        float diff = velocityHistory[idx].velocity - mean;
        variance += diff * diff;
    }
    variance /= count;

    return sqrt(variance) < Config::VELOCITY_NOISE_THRESHOLD;
}

// Determines if the motor rotation is stopping based on velocity trend and magnitude
bool isRotationStopping() {
    float current_vel_magnitude = fabsf(smoothed_velocity);
    float trend_magnitude = fabsf(velocity_trend);
    float current_trend_sign = (velocity_trend > 0) ? 1.0f : -1.0f;

    // Rotation is considered stopping if:
    // 1. Velocity magnitude is relatively small (below consistency threshold)
    // 2. Velocity trend is negative (slowing down) AND its magnitude is significant
    // 3. The current velocity is consistent (not just noisy fluctuations)
    return (current_vel_magnitude < Config::MIN_VELOCITY_CONSISTENCY &&
            current_trend_sign < 0 &&
            trend_magnitude > Config::VELOCITY_TREND_THRESHOLD &&
            isVelocityConsistent(smoothed_velocity));
}

// Updates revolution count based on normalized angle difference
void updateRevolutionCount(float angle_diff) {
    // A full revolution is 2*PI. If angle_diff crosses PI, it indicates a wrap-around.
    // The threshold ensures we count a revolution only when it's clearly passed the half-point of a circle.
    const float REVOLUTION_WRAP_THRESHOLD = PI * Config::REVOLUTION_THRESHOLD_FACTOR;

    if (angle_diff > REVOLUTION_WRAP_THRESHOLD) {
        revolution_count--; // Angle unwrapped from positive PI to negative PI
    } else if (angle_diff < -REVOLUTION_WRAP_THRESHOLD) {
        revolution_count++; // Angle unwrapped from negative PI to positive PI
    }
}

void applyEKF(float measured_angle) {
    unsigned long now = micros();
    // Calculate delta time in seconds
    float dt = (float)(now - prev_ekf_micros) / 1e6f;
    prev_ekf_micros = now;

    // --- Prediction Step ---
    float angle = x_ekf(0);
    float velocity = x_ekf(1);

    // Predict next state: angle = angle + velocity * dt, velocity remains constant
    float predicted_angle = angle + velocity * dt; // Angle not normalized yet for prediction
    float predicted_velocity = velocity;

    BLA::Matrix<N_STATES_EKF, 1> x_pred = {predicted_angle, predicted_velocity};

    // Calculate Jacobian of state transition (A)
    // A = d(f)/dx = [1 dt; 0 1]
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> A = {{1.0f, dt}, {0.0f, 1.0f}};

    // Predict error covariance: P_pred = A * P_ekf * A_transpose + Q_ekf
    P_ekf = A * P_ekf * BLA::Transpose(A) + Q_ekf;

    // --- Update Step ---
    // Measurement function H maps state to measurement: H = [1 0] (only angle is measured directly)
    // Innovation (measurement residual): y = z - H * x_pred
    // Normalize the innovation (angle difference)
    float y = normalizeAngle(measured_angle - predicted_angle);

    // Innovation covariance: S = H * P_pred * H_transpose + R_ekf
    // Since H = [1 0], H * P_pred * H_transpose simplifies to P_pred(0,0)
    BLA::Matrix<1, 1> S = R_ekf + P_ekf(0, 0);

    // Kalman Gain: K = P_pred * H_transpose * S_inverse
    // K = [P_pred(0,0)/S(0,0); P_pred(1,0)/S(0,0)]
    BLA::Matrix<N_STATES_EKF, 1> K = {P_ekf(0, 0) / S(0, 0), P_ekf(1, 0) / S(0, 0)};

    // Update state estimate: x_ekf = x_pred + K * y
    x_ekf = x_pred + K * y;
    x_ekf(0) = normalizeAngle(x_ekf(0)); // Normalize angle after update

    // Update error covariance: P_ekf = (I - K * H) * P_pred
    // (I - K * H) simplifies to I_matrix_ekf with K(0,0) and K(1,0) subtracted from first column
    BLA::Matrix<N_STATES_EKF, N_STATES_EKF> I_KH;
    I_KH = I_matrix_ekf;
    I_KH(0, 0) -= K(0, 0);
    I_KH(1, 0) -= K(1, 0); // K(1,0) is equivalent to P_pred(1,0)/S(0,0)
    P_ekf = I_KH * P_ekf;

    // --- Update History and Tracking Variables ---
    updateVelocityHistory(x_ekf(1), now);
    updateAngleHistory(x_ekf(0), now); // Keep track of estimated angle
    smoothed_velocity = getSmoothedVelocity();
    velocity_trend = getVelocityTrend();
    is_stopping = isRotationStopping();

    // --- Enhanced EKF Locking Logic ---
    float current_vel_magnitude = fabsf(x_ekf(1)); // Use EKF's estimated velocity directly for 'is_moving'
    bool is_moving_threshold = current_vel_magnitude > Config::VELOCITY_STEP_THRESHOLD;
    bool current_is_consistent = isVelocityConsistent(x_ekf(1)); // Check consistency of current EKF velocity

    // Direction consistency check
    if (is_moving_threshold && current_is_consistent) {
        float current_direction_sign = (x_ekf(1) > 0) ? 1.0f : -1.0f;
        float last_stable_direction_sign = (last_stable_velocity > 0) ? 1.0f : -1.0f;

        if (fabsf(last_stable_velocity) > Config::VELOCITY_STEP_THRESHOLD && current_direction_sign == last_stable_direction_sign) {
            consistent_direction_count++;
        } else {
            consistent_direction_count = 0;
            last_direction_change_micros = now;
        }
        last_stable_velocity = x_ekf(1); // Update last stable velocity
    } else {
        consistent_direction_count = 0; // Reset if not moving or inconsistent
    }

    // --- EKF Lock/Unlock State Management ---
    if (!ekfLocked) {
        // Attempt to lock if conditions met
        if (is_moving_threshold && current_is_consistent && consistent_direction_count >= Config::CONSISTENT_DIR_COUNT_THRESHOLD) {
            ekfLocked = true;
            Serial.println(F("EKF Locked: Consistent movement detected."));
        }
    } else { // EKF is currently locked
        // Update unwrapped angle when locked
        float angle_diff = normalizeAngle(x_ekf(0) - unwrapped_electrical_angle);
        unwrapped_electrical_angle += angle_diff;
        updateRevolutionCount(angle_diff);
        lastNonZeroVelMicros = now; // Reset timeout as long as locked

        // Check for unlock conditions
        bool should_unlock = false;
        if (is_stopping) { // Actively stopping detected
            should_unlock = true;
            Serial.println(F("EKF Unlocked: Rotation stopping detected."));
        } else if (!is_moving_threshold && (now - lastNonZeroVelMicros > (unsigned long)Config::ZERO_VEL_LOCK_TIMEOUT_MS * 1000UL)) {
            // Prolonged zero velocity and not moving
            should_unlock = true;
            Serial.println(F("EKF Unlocked: Prolonged zero velocity timeout."));
        } else if (!current_is_consistent && current_vel_magnitude < Config::MIN_VELOCITY_CONSISTENCY) {
            // Inconsistent low velocity
            should_unlock = true;
            Serial.println(F("EKF Unlocked: Inconsistent low velocity."));
        }

        if (should_unlock) {
            ekfLocked = false;
            consistent_direction_count = 0; // Reset for next lock
            // Optionally, reset EKF velocity to zero or a small value upon unlock
            // x_ekf(1) = 0.0f;
            // P_ekf(1,1) = 1.0f; // Increase velocity uncertainty
        }
    }
    prev_is_moving = is_moving_threshold; // Track previous state for transitions
}

// Normalizes an angle to the range (-PI, PI]
float normalizeAngle(float angle) {
    while (angle <= -PI) angle += Config::TWO_PI_F;
    while (angle > PI) angle -= Config::TWO_PI_F;
    return angle;
}

// Safely reads analog pin, clamping values to 0-1023
int checkedAnalogRead(int pin) {
    int val = analogRead(pin);
    // These checks are mostly for defensive programming; analogRead typically stays within bounds.
    if (val < 0) return 0;
    if (val > 1023) return 1023;
    return val;
}

// Logs various state and estimation data to Serial
void logData() {
    // Using F() macro for strings to save RAM
    Serial.print(F("Ang: ")); Serial.print(x_ekf(0), 3);
    Serial.print(F("\tVel: ")); Serial.print(x_ekf(1), 4); // EKF estimated velocity
    Serial.print(F("\tSmVel: ")); Serial.print(smoothed_velocity, 4); // Smoothed history velocity
    Serial.print(F("\tTrend: ")); Serial.print(velocity_trend, 4);
    Serial.print(F("\tSteps: ")); Serial.print(estimated_mechanical_steps);
    Serial.print(F("\tRev: ")); Serial.print(revolution_count);
    Serial.print(F("\tLock: ")); Serial.print(ekfLocked ? F("Y") : F("N"));
    Serial.print(F("\tStop: ")); Serial.println(is_stopping ? F("Y") : F("N"));
    // Serial.print(F("\tDirCnt: ")); Serial.print(consistent_direction_count);
    // Serial.print(F("\tP(0,0): ")); Serial.print(P_ekf(0,0), 6); // Debug: EKF Angle Covariance
    // Serial.print(F("\tP(1,1): ")); Serial.print(P_ekf(1,1), 6); // Debug: EKF Velocity Covariance
}
