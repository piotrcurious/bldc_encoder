#include <Arduino.h>
#include <BasicLinearAlgebra.h> // Include the BLA library

// Uncomment this line to enable DIAGNOSTIC MODE:
// In this mode, the code will continuously read and print RAW ADC values for A0, A1, A2
// without any motor control logic or auto-tuning, allowing you to directly
// troubleshoot your hardware setup with a multimeter.
//#define DIAGNOSTIC_MODE

const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int COMMON_PIN = 9; // Pin 9 is controlled by Timer1 (along with pin 10)

volatile long encoderPosition = 0;
int previousHallState = 0;
unsigned long lastStateChangeMillis = 0;

const unsigned long DEBOUNCE_DELAY_MS = 1;
const int NUM_SAMPLES = 16;
const float VOLTAGE_THRESHOLD = 0.08f; // Reduced threshold for better low RPM sensitivity

// Excitation parameters
int EXCITATION_PWM_VALUE = 1; // THIS WILL BE AUTO-TUNED. Initial guess for PWM strength (0-255).
const int EXCITATION_PULSE_WIDTH_US = 1000; // Duration of the initial PWM pulse on COMMON_PIN.
const int FLOATING_SETTLING_TIME_US = 10;     // Crucial: Time to float after excitation pulse before sensing.

// Global variables for plotting (current averaged ADC values)
float g_valA = 0.0f;
float g_valB = 0.0f;
float g_valC = 0.0f;

// --- Enhanced Kalman Filter Variables ---
#define N_STATES 9 // Extended states: [V_A, V_B, V_C, dV_A/dt, dV_B/dt, dV_C/dt, d²V_A/dt², d²V_B/dt², d²V_C/dt²]
#define HISTORY_SIZE 32 // Circular buffer size for trend analysis
#define MIN_DIRECTION_CONFIDENCE 0.7f // Minimum confidence for direction change acceptance

// Enhanced state vector with velocity and acceleration components
BLA::Matrix<N_STATES, 1> x_hat;

// History buffers for trend analysis
struct HistoryBuffer {
    float values[HISTORY_SIZE];
    unsigned long timestamps[HISTORY_SIZE];
    int head;
    int count;
    
    void init() {
        head = 0;
        count = 0;
        for(int i = 0; i < HISTORY_SIZE; i++) {
            values[i] = 0.0f;
            timestamps[i] = 0;
        }
    }
    
    void push(float value, unsigned long timestamp) {
        values[head] = value;
        timestamps[head] = timestamp;
        head = (head + 1) % HISTORY_SIZE;
        if(count < HISTORY_SIZE) count++;
    }
    
    float getAverage() {
        if(count == 0) return 0.0f;
        float sum = 0.0f;
        for(int i = 0; i < count; i++) {
            sum += values[i];
        }
        return sum / count;
    }
    
    float getTrend() {
        if(count < 3) return 0.0f;
        
        // Calculate linear regression slope
        float sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        unsigned long base_time = timestamps[(head - count + HISTORY_SIZE) % HISTORY_SIZE];
        
        for(int i = 0; i < count; i++) {
            int idx = (head - count + i + HISTORY_SIZE) % HISTORY_SIZE;
            float x = (timestamps[idx] - base_time) / 1000.0f; // Convert to seconds
            float y = values[idx];
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }
        
        float denominator = count * sum_x2 - sum_x * sum_x;
        if(abs(denominator) < 1e-6) return 0.0f;
        
        return (count * sum_xy - sum_x * sum_y) / denominator;
    }
};

HistoryBuffer historyA, historyB, historyC;
HistoryBuffer stateHistory; // For tracking Hall state transitions

// Direction confidence tracking
float directionConfidence = 1.0f;
int consecutiveDirectionChanges = 0;
unsigned long lastDirectionChangeTime = 0;

// Enhanced State Transition Matrix with time-based dynamics
float dt = 0.001f; // Time step in seconds (will be updated dynamically)
BLA::Matrix<N_STATES, N_STATES> F;

// Observation Matrix - we only directly observe the voltage states
BLA::Matrix<3, N_STATES> H;

// Enhanced Process Noise Covariance Matrix (Q) - adaptive
BLA::Matrix<N_STATES, N_STATES> Q;

// Measurement Noise Covariance Matrix (R) - adaptive
BLA::Matrix<3, 3> R;

// Estimate Covariance Matrix (P)
BLA::Matrix<N_STATES, N_STATES> P;

// Identity matrices
BLA::Matrix<N_STATES, N_STATES> I_matrix;
BLA::Matrix<3, 3> I3_matrix;

// Measurement vector
BLA::Matrix<3, 1> z_k;

// Initialize identity matrices
void initialize_matrices() {
    // Initialize I_matrix (9x9 identity)
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            I_matrix(i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
    
    // Initialize I3_matrix (3x3 identity)
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            I3_matrix(i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
    
    // Initialize H matrix (3x9) - only observe voltage states
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            H(i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
}

// Update state transition matrix based on current time step
void updateStateTransitionMatrix(float delta_t) {
    dt = delta_t;
    
    // Reset F to identity
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            F(i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
    
    // Add velocity integration: position += velocity * dt
    for(int i = 0; i < 3; i++) {
        F(i, i + 3) = dt;
    }
    
    // Add acceleration integration: velocity += acceleration * dt
    for(int i = 3; i < 6; i++) {
        F(i, i + 3) = dt;
    }
    
    // Add second-order terms: position += 0.5 * acceleration * dt²
    for(int i = 0; i < 3; i++) {
        F(i, i + 6) = 0.5f * dt * dt;
    }
}

// Adaptive noise matrix calculation based on motion characteristics
void updateNoiseMatrices() {
    // Calculate motion intensity from velocity states
    float velocity_magnitude = sqrt(x_hat(3)*x_hat(3) + x_hat(4)*x_hat(4) + x_hat(5)*x_hat(5));
    float acceleration_magnitude = sqrt(x_hat(6)*x_hat(6) + x_hat(7)*x_hat(7) + x_hat(8)*x_hat(8));
    
    // Adaptive process noise - higher during rapid changes
    float base_q = 0.1f;
    float velocity_factor = 1.0f + velocity_magnitude * 0.1f;
    float acceleration_factor = 1.0f + acceleration_magnitude * 0.05f;
    float adaptive_q = base_q * velocity_factor * acceleration_factor;
    
    // Update Q matrix
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            if(i == j) {
                if(i < 3) Q(i, j) = adaptive_q; // Position states
                else if(i < 6) Q(i, j) = adaptive_q * 10.0f; // Velocity states
                else Q(i, j) = adaptive_q * 100.0f; // Acceleration states
            } else {
                Q(i, j) = 0.0f;
            }
        }
    }
    
    // Adaptive measurement noise - lower during stable conditions
    float stability_factor = 1.0f / (1.0f + velocity_magnitude * 0.1f);
    float adaptive_r = 50.0f * stability_factor;
    
    // Update R matrix
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            R(i, j) = (i == j) ? adaptive_r : 0.0f;
        }
    }
}

// Calculate direction confidence based on historical data and current state
float calculateDirectionConfidence(int currentState, int previousState) {
    if(currentState == 0 || previousState == 0) return 0.0f;
    
    // Get trends from history buffers
    float trendA = historyA.getTrend();
    float trendB = historyB.getTrend();
    float trendC = historyC.getTrend();
    
    // Calculate expected trends for the state transition
    float expectedTrendMagnitude = 0.0f;
    bool validTransition = false;
    
    // Check if this is a valid state transition
    int expectedNextStateCW[] = {0, 2, 3, 4, 5, 6, 1};
    int expectedNextStateCCW[] = {0, 6, 1, 2, 3, 4, 5};
    
    if(currentState == expectedNextStateCW[previousState]) {
        validTransition = true;
        expectedTrendMagnitude = 1.0f; // CW rotation
    } else if(currentState == expectedNextStateCCW[previousState]) {
        validTransition = true;
        expectedTrendMagnitude = -1.0f; // CCW rotation
    }
    
    if(!validTransition) return 0.0f;
    
    // Calculate consistency of trends with expected direction
    float trendConsistency = abs(trendA + trendB + trendC) / 3.0f;
    
    // Factor in velocity consistency from Kalman filter
    float velocityConsistency = sqrt(x_hat(3)*x_hat(3) + x_hat(4)*x_hat(4) + x_hat(5)*x_hat(5));
    
    // Combine factors
    float confidence = min(1.0f, trendConsistency * 0.3f + velocityConsistency * 0.7f);
    
    return confidence;
}

/**
 * @brief Applies an enhanced multivariate Kalman filter with velocity and acceleration states.
 */
void applyEnhancedKalmanFilter(int raw_measurements[], unsigned long timestamp) {
    static unsigned long lastTimestamp = 0;
    
    // Calculate time step
    if(lastTimestamp > 0) {
        dt = (timestamp - lastTimestamp) / 1000000.0f; // Convert microseconds to seconds
        dt = max(0.0001f, min(dt, 0.1f)); // Clamp dt to reasonable range
    }
    lastTimestamp = timestamp;
    
    // Update state transition matrix
    updateStateTransitionMatrix(dt);
    
    // Update adaptive noise matrices
    updateNoiseMatrices();
    
    // Convert raw measurements to BLA matrix
    z_k(0) = (float)raw_measurements[0];
    z_k(1) = (float)raw_measurements[1];
    z_k(2) = (float)raw_measurements[2];
    
    // --- Prediction Step ---
    // Predicted State: x_hat_k_pred = F * x_hat_k_minus_1
    BLA::Matrix<N_STATES, 1> x_hat_pred = F * x_hat;
    
    // Predicted Covariance: P_k_pred = F * P * F' + Q
    BLA::Matrix<N_STATES, N_STATES> F_transposed;
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            F_transposed(i, j) = F(j, i);
        }
    }
    BLA::Matrix<N_STATES, N_STATES> P_pred = F * P * F_transposed + Q;
    
    // --- Update Step ---
    // Innovation: y_k = z_k - H * x_hat_k_pred
    BLA::Matrix<3, 1> y_k = z_k - (H * x_hat_pred);
    
    // Innovation Covariance: S_k = H * P_k_pred * H' + R
    BLA::Matrix<3, N_STATES> H_transposed;
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            H_transposed(i, j) = H(j, i);
        }
    }
    BLA::Matrix<3, 3> S_k = H * P_pred * H_transposed + R;
    
    // Kalman Gain: K_k = P_k_pred * H' * S_k_inv
    BLA::Matrix<3, 3> S_k_inv = BLA::Inverse(S_k);
    
    // Check for inversion success
    if (isnan(S_k_inv(0,0))) {
        return; // Skip update if S_k is singular
    }
    
    BLA::Matrix<N_STATES, 3> K_k = P_pred * H_transposed * S_k_inv;
    
    // Updated State: x_hat_k = x_hat_k_pred + K_k * y_k
    x_hat = x_hat_pred + (K_k * y_k);
    
    // Updated Covariance: P_k = (I - K_k * H) * P_k_pred
    P = (I_matrix - K_k * H) * P_pred;
    
    // Update history buffers
    historyA.push(x_hat(0), timestamp);
    historyB.push(x_hat(1), timestamp);
    historyC.push(x_hat(2), timestamp);
}

// Auto-tuning parameters
const int TARGET_QUIESCENT_AVG_ADC = 80;
const int PWM_ADJUST_STEP = 1;
const int PWM_CALIBRATION_ITERATIONS = 2000;
const int PWM_CALIBRATION_TOLERANCE = 1;

/**
 * @brief Performs a single analog read.
 */
int reliableAnalogRead(int pin) {
    return analogRead(pin);
}

/**
 * @brief Performs a single cycle of PWM excitation, floating, and sensing.
 */
void executeSensingCycle() {
    unsigned long startTime = micros();
    
    // 1. Apply PWM excitation to COMMON_PIN
    analogWrite(COMMON_PIN, EXCITATION_PWM_VALUE);
    delayMicroseconds(EXCITATION_PULSE_WIDTH_US);

    // 2. Set COMMON_PIN to INPUT for BEMF sensing
    pinMode(COMMON_PIN, INPUT);
    delayMicroseconds(FLOATING_SETTLING_TIME_US);

    // 3. Aggressive oversampling
    long sumA = 0, sumB = 0, sumC = 0;
    for (int i = 0; i < NUM_SAMPLES; i++) {
        sumA += reliableAnalogRead(PHASE_A_PIN);
        sumB += reliableAnalogRead(PHASE_B_PIN);
        sumC += reliableAnalogRead(PHASE_C_PIN);
    }
    
    // Get raw averaged values
    int raw_vals[3];
    raw_vals[0] = sumA / NUM_SAMPLES;
    raw_vals[1] = sumB / NUM_SAMPLES;
    raw_vals[2] = sumC / NUM_SAMPLES;

    // Apply Enhanced Kalman Filter
    applyEnhancedKalmanFilter(raw_vals, micros());

    // Update global variables with filtered estimates
    g_valA = x_hat(0); 
    g_valB = x_hat(1);
    g_valC = x_hat(2);

    // 4. Re-establish default state
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);
}

/**
 * @brief Enhanced Hall state detection with confidence tracking.
 */
int getHallStateWithConfidence(float* confidence) {
    executeSensingCycle();
    
    // Determine state with enhanced sensitivity
    int detectedState = 0;
    float maxDifference = 0.0f;
    
    // Check all possible states and find the one with maximum confidence
    float states_confidence[7] = {0}; // Index 0 unused, states 1-6
    
    if (g_valA > (g_valB + VOLTAGE_THRESHOLD) && g_valB > (g_valC + VOLTAGE_THRESHOLD)) {
        states_confidence[1] = min(g_valA - g_valB, g_valB - g_valC);
    }
    if (g_valA > (g_valC + VOLTAGE_THRESHOLD) && g_valC > (g_valB + VOLTAGE_THRESHOLD)) {
        states_confidence[2] = min(g_valA - g_valC, g_valC - g_valB);
    }
    if (g_valB > (g_valA + VOLTAGE_THRESHOLD) && g_valA > (g_valC + VOLTAGE_THRESHOLD)) {
        states_confidence[3] = min(g_valB - g_valA, g_valA - g_valC);
    }
    if (g_valB > (g_valC + VOLTAGE_THRESHOLD) && g_valC > (g_valA + VOLTAGE_THRESHOLD)) {
        states_confidence[4] = min(g_valB - g_valC, g_valC - g_valA);
    }
    if (g_valC > (g_valA + VOLTAGE_THRESHOLD) && g_valA > (g_valB + VOLTAGE_THRESHOLD)) {
        states_confidence[5] = min(g_valC - g_valA, g_valA - g_valB);
    }
    if (g_valC > (g_valB + VOLTAGE_THRESHOLD) && g_valB > (g_valA + VOLTAGE_THRESHOLD)) {
        states_confidence[6] = min(g_valC - g_valB, g_valB - g_valA);
    }
    
    // Find state with highest confidence
    for(int i = 1; i <= 6; i++) {
        if(states_confidence[i] > maxDifference) {
            maxDifference = states_confidence[i];
            detectedState = i;
        }
    }
    
    *confidence = maxDifference / 10.0f; // Normalize confidence
    return detectedState;
}

/**
 * @brief Standard Hall state detection (backward compatibility).
 */
int getHallState() {
    float confidence;
    return getHallStateWithConfidence(&confidence);
}

/**
 * @brief Auto-tunes the EXCITATION_PWM_VALUE.
 */
void autoTuneExcitationPwm() {
    Serial.println("\n--- Starting Enhanced Auto-Tune ---");
    Serial.println("  Initializing enhanced Kalman filter...");
    
    // Initialize history buffers
    historyA.init();
    historyB.init();
    historyC.init();
    stateHistory.init();
    
    // Initialize enhanced Kalman filter
    for(int i = 0; i < N_STATES; ++i) {
        x_hat(i) = (i < 3) ? (float)TARGET_QUIESCENT_AVG_ADC : 0.0f;
        for(int j = 0; j < N_STATES; ++j) {
            P(i, j) = (i == j) ? 1000.0f : 0.0f;
        }
    }

    for (int i = 0; i < PWM_CALIBRATION_ITERATIONS; i++) {
        executeSensingCycle();

        int current_avg_adc = (int)((g_valA + g_valB + g_valC) / 3.0f);

        Serial.print("Enhanced Tune Step "); Serial.print(i + 1);
        Serial.print("\tPWM: "); Serial.print(EXCITATION_PWM_VALUE);
        Serial.print("\tAvg ADC: "); Serial.print(current_avg_adc);
        Serial.print("\tVelocity: "); Serial.print(sqrt(x_hat(3)*x_hat(3) + x_hat(4)*x_hat(4) + x_hat(5)*x_hat(5)));
        Serial.print("\tPhases (A/B/C): "); Serial.print(g_valA);
        Serial.print("/"); Serial.print(g_valB);
        Serial.print("/"); Serial.println(g_valC);

        if (abs(current_avg_adc - TARGET_QUIESCENT_AVG_ADC) <= PWM_CALIBRATION_TOLERANCE) {
            Serial.println("--- Enhanced auto-tune complete! ---");
            break;
        }

        if (current_avg_adc < TARGET_QUIESCENT_AVG_ADC) {
            EXCITATION_PWM_VALUE += PWM_ADJUST_STEP;
        } else {
            EXCITATION_PWM_VALUE -= PWM_ADJUST_STEP;
        }

        EXCITATION_PWM_VALUE = max(0, min(255, EXCITATION_PWM_VALUE));
        delay(50);
    }
    
    Serial.print("Final auto-tuned EXCITATION_PWM_VALUE: "); Serial.println(EXCITATION_PWM_VALUE);
    Serial.println("Enhanced Kalman filter initialized with history tracking.");
    Serial.println("--------------------------------------------------");
}

void setup() {
    Serial.begin(115200);
    Serial.println("Enhanced BLDC Encoder Initialized.");

    // Initialize matrices
    initialize_matrices();

    // ADC setup
    analogReference(INTERNAL);
    ADCSRA &= ~((1 << ADPS0) | (1 << ADPS1) | (1 << ADPS2));
    ADCSRA |= (1 << ADPS2) | (1 << ADPS0);
    ADCSRA |= (1 << ADEN);

    // Increase PWM frequency for pin 9
    TCCR1B = (TCCR1B & 0xF8) | 0x01;

    // Common pin setup
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);

#ifndef DIAGNOSTIC_MODE
    // Run enhanced auto-tuning
    autoTuneExcitationPwm();

    // Initialize state tracking
    previousHallState = getHallState();
    lastStateChangeMillis = millis();
    directionConfidence = 1.0f;
#endif
}

void loop() {
#ifdef DIAGNOSTIC_MODE
    int valA_raw = reliableAnalogRead(PHASE_A_PIN);
    int valB_raw = reliableAnalogRead(PHASE_B_PIN);
    int valC_raw = reliableAnalogRead(PHASE_C_PIN);

    Serial.print("RAW ADC (A/B/C): ");
    Serial.print(valA_raw); Serial.print("\t");
    Serial.print(valB_raw); Serial.print("\t");
    Serial.println(valC_raw);
    delay(100);
#else
    float stateConfidence;
    int currentHallState = getHallStateWithConfidence(&stateConfidence);

    if (currentHallState != previousHallState && currentHallState != 0) {
        // Calculate direction confidence
        float directionConf = calculateDirectionConfidence(currentHallState, previousHallState);
        
        // Only accept direction changes with sufficient confidence
        if (directionConf >= MIN_DIRECTION_CONFIDENCE || stateConfidence > 0.8f) {
            int expectedNextStateCW[] = {0, 2, 3, 4, 5, 6, 1};
            int expectedNextStateCCW[] = {0, 6, 1, 2, 3, 4, 5};

            if (currentHallState == expectedNextStateCW[previousHallState]) {
                encoderPosition++;
                directionConfidence = min(1.0f, directionConfidence + 0.1f);
            } else if (currentHallState == expectedNextStateCCW[previousHallState]) {
                encoderPosition--;
                directionConfidence = min(1.0f, directionConfidence + 0.1f);
            } else {
                // Unexpected state transition - reduce confidence
                directionConfidence = max(0.1f, directionConfidence - 0.2f);
            }
            
            previousHallState = currentHallState;
            lastStateChangeMillis = millis();
            
            // Update state history
            stateHistory.push((float)currentHallState, millis());
        }
    }

    // Check for serial input
    if (Serial.available()) {
        char command = Serial.read();
        if (command == 'T' || command == 't') {
            autoTuneExcitationPwm();
            encoderPosition = 0;
            previousHallState = getHallState();
            lastStateChangeMillis = millis();
            directionConfidence = 1.0f;
        }
    }

    // Enhanced output for Serial Plotter
    Serial.print(g_valA); Serial.print("\t");
    Serial.print(g_valB); Serial.print("\t");
    Serial.print(g_valC); Serial.print("\t");
    Serial.print(encoderPosition); Serial.print("\t");
    Serial.print(currentHallState); Serial.print("\t");
    Serial.print(directionConfidence); Serial.print("\t");
    Serial.print(sqrt(x_hat(3)*x_hat(3) + x_hat(4)*x_hat(4) + x_hat(5)*x_hat(5))); // Velocity magnitude
    Serial.print("\t");
    Serial.println(stateConfidence);

#endif
}
