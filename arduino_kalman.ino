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
// Changed VOLTAGE_THRESHOLD to float to leverage filtered accuracy
const float VOLTAGE_THRESHOLD = 0.1f;

// Excitation parameters
int EXCITATION_PWM_VALUE = 1; // THIS WILL BE AUTO-TUNED. Initial guess for PWM strength (0-255).
const int EXCITATION_PULSE_WIDTH_US = 1000; // Duration of the initial PWM pulse on COMMON_PIN.
const int FLOATING_SETTLING_TIME_US = 10;     // Crucial: Time to float after excitation pulse before sensing.
                                            // YOU MUST MANUALLY OPTIMIZE THIS WITH SERIAL PLOTTER.

// Global variables for plotting (current averaged ADC values)
// Changed g_valA, g_valB, g_valC to float to retain Kalman filter precision
float g_valA = 0.0f;
float g_valB = 0.0f;
float g_valC = 0.0f;

// --- Kalman Filter Variables ---
#define N_STATES 3 // Number of states (V_A, V_B, V_C)

// State vector [V_A, V_B, V_C]'
BLA::Matrix<N_STATES, 1> x_hat; // Will be initialized with a reasonable guess in setup()

// State Transition Matrix (F) - Identity matrix, assuming state doesn't change between samples
BLA::Matrix<N_STATES, N_STATES> F = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

// Observation Matrix (H) - Identity matrix, assuming measurements directly correspond to states
BLA::Matrix<N_STATES, N_STATES> H = {
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
};

// Process Noise Covariance Matrix (Q)
// Updated Q_DIAGONAL_VALUE to 1.0f as per user's working value
const float Q_DIAGONAL_VALUE = 1.0f; 
BLA::Matrix<N_STATES, N_STATES> Q = {
    Q_DIAGONAL_VALUE, 0, 0,
    0, Q_DIAGONAL_VALUE, 0,
    0, 0, Q_DIAGONAL_VALUE
};

// Measurement Noise Covariance Matrix (R)
// Updated R_DIAGONAL_VALUE to 100.0f as per user's working value
const float R_DIAGONAL_VALUE = 100.0f; 
BLA::Matrix<N_STATES, N_STATES> R = {
    R_DIAGONAL_VALUE, 0, 0,
    0, R_DIAGONAL_VALUE, 0,
    0, 0, R_DIAGONAL_VALUE
};

// Initial Estimate Covariance Matrix (P)
BLA::Matrix<N_STATES, N_STATES> P; // Initialized in setup() or autoTuneExcitationPwm()

// Manual initialization for Identity matrix (I_matrix)
BLA::Matrix<N_STATES, N_STATES> I_matrix;
// Initialize I_matrix manually since BLA::Identity() isn't recognized
void initialize_I_matrix() {
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            I_matrix(i, j) = (i == j) ? 1.0f : 0.0f;
        }
    }
}


// Measurement vector (z_k)
BLA::Matrix<N_STATES, 1> z_k;

/**
 * @brief Applies a multivariate Kalman filter to the phase measurements.
 * Updates the global x_hat (state estimate) and P (covariance) matrices.
 * @param raw_measurements An array of raw ADC readings [A, B, C].
 */
void applyMultivariateKalmanFilter(int raw_measurements[]) {
    // Convert raw measurements to BLA matrix
    z_k(0) = (float)raw_measurements[0];
    z_k(1) = (float)raw_measurements[1];
    z_k(2) = (float)raw_measurements[2];

    // --- Prediction Step ---
    // Predicted State: x_hat_k_pred = F * x_hat_k_minus_1
    // (Since F is identity, x_hat_k_pred remains x_hat. No explicit multiplication needed)

    // Predicted Covariance: P_k_pred = F * P * F' + Q
    // For F = Identity, this simplifies to P_k_pred = P + Q
    BLA::Matrix<N_STATES, N_STATES> P_pred = P + Q; // BLA simplifies this

    // --- Update Step ---
    // Innovation (Measurement Residual): y_k = z_k - H * x_hat_k_pred
    BLA::Matrix<N_STATES, 1> y_k = z_k - (H * x_hat);

    // Innovation Covariance: S_k = H * P_k_pred * H' + R
    // For H = Identity, this simplifies to S_k = P_pred + R
    BLA::Matrix<N_STATES, N_STATES> S_k = P_pred + R; // BLA simplifies this

    // Kalman Gain: K_k = P_k_pred * H' * S_k_inv
    // Attempting BLA::Inverse() and manual Transpose.
    BLA::Matrix<N_STATES, N_STATES> S_k_inv = BLA::Inverse(S_k);

    // Check for inversion success (BLA returns NaN if singular)
    if (isnan(S_k_inv(0,0))) { // Check if the first element is NaN (indicates inversion failure)
        // You might want to log this error or handle it differently in a real application
        return; // Skip update if S_k is singular
    }

    // Manual Transpose for H, as BLA::Transpose() is not recognized
    BLA::Matrix<N_STATES, N_STATES> H_transposed;
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            H_transposed(i, j) = H(j, i);
        }
    }

    BLA::Matrix<N_STATES, N_STATES> K_k = P_pred * H_transposed * S_k_inv;

    // Updated State: x_hat_k = x_hat_k_pred + K_k * y_k
    x_hat = x_hat + (K_k * y_k); // Use the predicted state as current x_hat for the addition

    // Updated Covariance: P_k = (I - K_k * H) * P_k_pred
    P = (I_matrix - K_k * H) * P_pred;
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
    int raw_vals[N_STATES];
    raw_vals[0] = sumA / NUM_SAMPLES;
    raw_vals[1] = sumB / NUM_SAMPLES;
    raw_vals[2] = sumC / NUM_SAMPLES;

    // Apply Multivariate Kalman Filter
    applyMultivariateKalmanFilter(raw_vals);

    // Update global variables with filtered estimates, now retaining float precision
    g_valA = x_hat(0); 
    g_valB = x_hat(1);
    g_valC = x_hat(2);

    // 4. Re-establish default state: COMMON_PIN as OUTPUT, LOW (stable ground sink)
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);
}

/**
 * @brief Reads analog voltages and determines Hall state based on collected data.
 * This function calls executeSensingCycle() to get the latest phase data.
 * @return Current Hall state (1-6) or 0 if ambiguous.
 */
int getHallState() {
    executeSensingCycle(); // Perform excitation and data collection

    // Determine state with sensitive hysteresis using FILTERED values (now float)
    if (g_valA > (g_valB + VOLTAGE_THRESHOLD) && g_valB > (g_valC + VOLTAGE_THRESHOLD)) return 1;
    if (g_valA > (g_valC + VOLTAGE_THRESHOLD) && g_valC > (g_valB + VOLTAGE_THRESHOLD)) return 2;
    if (g_valB > (g_valA + VOLTAGE_THRESHOLD) && g_valA > (g_valC + VOLTAGE_THRESHOLD)) return 3;
    if (g_valB > (g_valC + VOLTAGE_THRESHOLD) && g_valC > (g_valA + VOLTAGE_THRESHOLD)) return 4;
    if (g_valC > (g_valA + VOLTAGE_THRESHOLD) && g_valA > (g_valB + VOLTAGE_THRESHOLD)) return 5;
    if (g_valC > (g_valB + VOLTAGE_THRESHOLD) && g_valB > (g_valA + VOLTAGE_THRESHOLD)) return 6;

    return 0;
}

/**
 * @brief Auto-tunes the EXCITATION_PWM_VALUE to center the average phase voltage.
 * Call this function when the motor is stationary.
 */
void autoTuneExcitationPwm() {
    Serial.println("\n--- Starting Auto-Tune for EXCITATION_PWM_VALUE ---");
    Serial.println("  (IMPORTANT: Ensure motor is ABSOLUTELY STATIONARY!)");
    Serial.println("  Monitoring average phase voltage during tuning.");
    Serial.println("  Press 'T' in Serial Monitor to re-tune at any time.");
    Serial.println("--------------------------------------------------");

    // Re-initialize Kalman filter states and covariance matrix for a clean start
    // Manual initialization for P (as identity with scale)
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            P(i, j) = (i == j) ? 1000.0f : 0.0f; // Set diagonal elements for high uncertainty
        }
    }
    // Manual initialization for x_hat (zeros)
    for(int i = 0; i < N_STATES; ++i) {
        x_hat(i) = 0.0f;
    }

    for (int i = 0; i < PWM_CALIBRATION_ITERATIONS; i++) {
        executeSensingCycle(); // Get current sensor readings with current PWM value (now filtered)

        // Cast to int for average calculation to avoid float math over time
        int current_avg_adc = (int)((g_valA + g_valB + g_valC) / 3.0f);

        Serial.print("Tune Step "); Serial.print(i + 1);
        Serial.print("\tPWM: "); Serial.print(EXCITATION_PWM_VALUE);
        Serial.print("\tAvg ADC: "); Serial.print(current_avg_adc);
        Serial.print("\tPhases (A/B/C): "); Serial.print(g_valA);
        Serial.print("/"); Serial.print(g_valB);
        Serial.print("/"); Serial.println(g_valC);


        if (abs(current_avg_adc - TARGET_QUIESCENT_AVG_ADC) <= PWM_CALIBRATION_TOLERANCE) {
            Serial.println("--- Auto-tune complete! Optimal PWM found. ---");
            break;
        }

        if (current_avg_adc < TARGET_QUIESCENT_AVG_ADC) {
            EXCITATION_PWM_VALUE += PWM_ADJUST_STEP;
        } else {
            EXCITATION_PWM_VALUE -= PWM_ADJUST_STEP;
        }

        // Clamp PWM value within valid range (0-255)
        if (EXCITATION_PWM_VALUE < 0) EXCITATION_PWM_VALUE = 0;
        if (EXCITATION_PWM_VALUE > 255) EXCITATION_PWM_VALUE = 255;

        delay(50); // Small delay between calibration steps for stability
    }
    Serial.print("Final auto-tuned EXCITATION_PWM_VALUE: "); Serial.println(EXCITATION_PWM_VALUE);
    Serial.println("--------------------------------------------------");
    Serial.println("Resume operation. Rotate motor for data.");
    Serial.println(" "); // Blank line for readability
}


void setup() {
    Serial.begin(115200);
    Serial.println("BLDC Encoder Initialized.");

    // Call the function to initialize I_matrix once in setup
    initialize_I_matrix();

    // ADC setup (faster conversion, internal 1.1V reference)
    analogReference(INTERNAL);
    // Setting ADC prescaler to 32 for ~46.8 kHz ADC clock (16MHz/32 = 500kHz, 13 cycles/conversion, so ~38kHz conversion rate)
    ADCSRA &= ~((1 << ADPS0) | (1 << ADPS1) | (1 << ADPS2)); // Clear prescaler bits
    ADCSRA |= (1 << ADPS2) | (1 << ADPS0); // Set prescaler to 32 (16MHz / 32 = 500kHz ADC clock)
    ADCSRA |= (1 << ADEN); // Enable ADC

    // *** INCREASE PWM FREQUENCY FOR PIN 9 (Timer1) ***
    // Pins 9 and 10 are controlled by Timer1.
    // Default Arduino Uno/Nano PWM frequency is 490Hz (prescaler 64)
    // To increase it to 31.25 kHz (16MHz / (1 * 512)) for 8-bit PWM, we set prescaler to 1.
    TCCR1B = (TCCR1B & 0xF8) | 0x01; // Set CS12, CS11, CS10 to 001 for prescaler of 1

    // Common pin setup (default to output LOW, strong sink)
    pinMode(COMMON_PIN, OUTPUT);
    digitalWrite(COMMON_PIN, LOW);

#ifndef DIAGNOSTIC_MODE
    // Run auto-tuning at startup
    autoTuneExcitationPwm();

    // After tuning, initialize x_hat to a sensible starting point (e.g., target average)
    x_hat(0) = (float)TARGET_QUIESCENT_AVG_ADC;
    x_hat(1) = (float)TARGET_QUIESCENT_AVG_ADC;
    x_hat(2) = (float)TARGET_QUIESCENT_AVG_ADC;
    
    // Reset P to reflect that x_hat is now based on a measurement, so initial uncertainty is lower.
    // Manual initialization for P (as identity with scale)
    for(int i = 0; i < N_STATES; ++i) {
        for(int j = 0; j < N_STATES; ++j) {
            P(i, j) = (i == j) ? 100.0f : 0.0f; // Lower initial uncertainty after initial guess
        }
    }

    previousHallState = getHallState(); // Get Hall state based on this new x_hat
    lastStateChangeMillis = millis();
#endif
}

void loop() {
#ifdef DIAGNOSTIC_MODE
    // Diagnostic Mode: Read and print raw ADC values for debugging hardware
    // This will *not* use the excitation or Kalman filter.
    int valA_raw = reliableAnalogRead(PHASE_A_PIN);
    int valB_raw = reliableAnalogRead(PHASE_B_PIN);
    int valC_raw = reliableAnalogRead(PHASE_C_PIN);

    Serial.print("RAW ADC (A/B/C): ");
    Serial.print(valA_raw); Serial.print("\t");
    Serial.print(valB_raw); Serial.print("\t");
    Serial.println(valC_raw);
    delay(100); // Slow down for easier observation in diagnostic mode
#else
    int currentHallState = getHallState();

    if (currentHallState != previousHallState && currentHallState != 0) {
        int expectedNextStateCW[] = {0, 2, 3, 4, 5, 6, 1};
        int expectedNextStateCCW[] = {0, 6, 1, 2, 3, 4, 5};

        if (currentHallState == expectedNextStateCW[previousHallState]) {
            encoderPosition++;
        } else if (currentHallState == expectedNextStateCCW[previousHallState]) {
            encoderPosition--;
        }
        previousHallState = currentHallState;
        lastStateChangeMillis = millis();
    }

    // Check for serial input to re-tune
    if (Serial.available()) {
        char command = Serial.read();
        if (command == 'T' || command == 't') {
            autoTuneExcitationPwm(); // Re-run tuning if 'T' is sent
            // Reset position and state after re-tuning for clarity
            encoderPosition = 0;
            // Re-initialize x_hat and P after tuning, similar to setup
            // x_hat is updated inside executeSensingCycle already based on initial average
            x_hat(0) = (float)TARGET_QUIESCENT_AVG_ADC;
            x_hat(1) = (float)TARGET_QUIESCENT_AVG_ADC;
            x_hat(2) = (float)TARGET_QUIESCENT_AVG_ADC;
            // Manual initialization for P (as identity with scale)
            for(int i = 0; i < N_STATES; ++i) {
                for(int j = 0; j < N_STATES; ++j) {
                    P(i, j) = (i == j) ? 100.0f : 0.0f;
                }
            }
            previousHallState = getHallState();
            lastStateChangeMillis = millis();
        }
    }

    // Output for Serial Plotter: valA, valB, valC, encoderPosition, currentHallState
    // These are now FILTERED values (floats)
    Serial.print(g_valA);
    Serial.print("\t");
    Serial.print(g_valB);
    Serial.print("\t");
    Serial.print(g_valC);
    Serial.print("\t");
    Serial.print(encoderPosition);
    Serial.print("\t");
    Serial.println(currentHallState);

//    delay(1); // Small delay if needed for stable serial output, but remove for performance
#endif
}
