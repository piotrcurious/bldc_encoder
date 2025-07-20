#include <Arduino.h>
#include <BasicLinearAlgebra.h>
#include <limits.h>

// --- Configuration ---
namespace Config {
    // Hardware Pins
    constexpr int PA = A0, PB = A1, PC = A2, CP = 9;

    // Motor & Angle Constants
    constexpr float ESR = 6.0f; // ELECTRICAL_STEPS_PER_REV
    const float PF = 3.14159265358979323846f; // PI_F
    constexpr float TPF = 2.0f * PI; // TWO_PI_F
    constexpr float SQ3 = 1.7320508f; // SQRT_3

    // ADC Sampling
    constexpr int DNS = 4; // DEFAULT_NUM_SAMPLES
    constexpr int MIN_S = 1; // MIN_NUM_SAMPLES
    constexpr int MAX_S = 16; // MAX_NUM_SAMPLES

    // Magnitude Thresholds for Dynamic Oversampling
    constexpr float LOW_MAGNITUDE_THRESHOLD = 5.0f; // Below this, use MAX_S samples
    constexpr float HIGH_MAGNITUDE_THRESHOLD = 7.0f; // Above this, use MIN_S samples

    // Dithering Configuration
    constexpr int DITHER_PULSE_WIDTH_US = 32; // Duration of the dither pulse in microseconds
    constexpr int DITHER_AMPLITUDE_PIN_STATE = 32; // Or LOW, or an analogWrite value for more subtle dither

    // Excitation & Auto-Tuning
    constexpr int EPWU = 31, FSTU = 1;
    constexpr int TQA = 32, PAS = 1, PCI = 2000, PCT = 1;

    // Kalman Filter (for ADC smoothing)
    constexpr float QDK1_BASE = 1.0f; // Base process noise
    constexpr float RDK1_BASE = 20.1f; // Base measurement noise
    constexpr float MATRIX_SINGULARITY_THRESHOLD = 1e-3f;

    // --- Core Estimation & Tracking Parameters (REVISED) ---
    constexpr float MIN_SIGNAL_MAGNITUDE = 0.08f; // Minimum vector magnitude for a reliable signal (in ADC units).

    // Velocity Smoothing (Exponential Moving Average)
    constexpr float VELOCITY_ALPHA = 0.2f; // EMA alpha for velocity (0.0 to 1.0, higher means less smoothing)
    constexpr float VELOCITY_ENGAGE_THRESHOLD = 0.05f; // Velocity magnitude to start increasing confidence
    constexpr float VELOCITY_DISENGAGE_THRESHOLD = 0.2f; // Velocity magnitude to start decreasing confidence (lower than engage)

    // Confidence-based State Transitions
    constexpr float CONFIDENCE_GAIN_RATE = 0.005f; // How quickly confidence increases per reliable frame
    constexpr float CONFIDENCE_DECAY_RATE = 0.02f; // How quickly confidence decreases per unreliable frame
    constexpr float ENGAGE_CONFIDENCE_THRESHOLD = 0.4f; // Confidence level to transition to ENGAGING
    constexpr float LOCK_CONFIDENCE_THRESHOLD = 0.8f;    // Confidence level to transition to LOCKED
    constexpr float UNLOCK_CONFIDENCE_THRESHOLD = 0.2f; // Confidence level to transition back to UNLOCKED

    // --- Unbalance Detection Thresholds --- (NEW)
    constexpr float ZERO_COMPONENT_UNBALANCE_THRESHOLD = 5.0f; // Threshold for zero component to indicate unbalance (in ADC units)
    constexpr float DQ_UNBALANCE_THRESHOLD = 5.0f; // Threshold for d/q components in specific scenarios (e.g., stationary)
    constexpr float KALMAN_RESIDUAL_THRESHOLD = 10.0f; // Sum of absolute residuals to indicate high measurement deviation
    constexpr float PHASE_MAGNITUDE_IMBALANCE_RATIO = 0.1f; // e.g., 0.1 means 10% deviation allowed between max/min phase magnitude

}

// --- Global State Variables ---

// Kalman Filter for ADC reading stabilization
#define NKF1 3
BLA::Matrix<NKF1, 1> xk1;
BLA::Matrix<NKF1, NKF1> Pk1, Qk1, Rk1, Ik1;

// Excitation PWM value (tuned at runtime)
int EPW = 1;

// Core State Machine
enum TrackingState { UNLOCKED, ENGAGING, LOCKED };
TrackingState currentState = UNLOCKED;
float trackingConfidence = 0.0f; // NEW: Continuous confidence level (0.0 to 1.0)

// Angle & Velocity State
float lastRawElectricalAngle = 0.0f; // Stores the raw angle from the *previous* reading, used by updateUnwrappedMeasuredAngle.
float unwrappedMeasuredAngle = 0.0f; // The continuously tracked electrical angle.
float estimatedAngle = 0.0f; // Output angle, usually unwrappedMeasuredAngle when locked.
float estimatedVelocity = 0.0f; // Output velocity (EMA smoothed).

// Store the unwrapped angle when we lose lock
float last_known_unwrapped_angle = 0.0f;
bool has_been_locked_before = false; // Flag to indicate if we've ever achieved a lock

// Timestamps and values for velocity calculation
unsigned long pVTs = 0; // Previous Value Timestamp (used for logging, less critical for core vel calc)
unsigned long last_velocity_update_timestamp = 0; // Timestamp of the last *successful* velocity calculation.
float last_unwrapped_angle_for_velocity_calc = 0.0f; // The unwrapped angle used in the *previous* velocity calculation.

// Mechanical Position State
long rC = 0;  // revolution_count
long eMS = 0; // estimated_mechanical_steps

// --- Structs ---
struct PhaseAnalysisResult {
    float angle = 0.0f;
    float magnitude = 0.0f;
    bool is_reliable = false;
    // New members for DQ0 transform results
    float d_component = 0.0f;
    float q_component = 0.0f;
    float zero_component = 0.0f;
    // New members for unbalance detection
    float kalman_residual_sum_abs = 0.0f; // Sum of absolute residuals from Kalman filter
    bool is_unbalanced = false; // Flag indicating detected unbalance
};

// --- Function Prototypes ---
void initM();
void aTEP();
void eSC(float current_signal_magnitude); // Modified: Added parameter
PhaseAnalysisResult analyzePhases();
void updateUnwrappedMeasuredAngle(float currentRawAngle);
void lD(const PhaseAnalysisResult& result);
int cAR(int); // Modified for dithering
void aMKF(const BLA::Matrix<NKF1, 1>& z);
int gDNS(float currentMagnitude); // Modified: Changed parameter
unsigned long safeMicrosDiff(unsigned long current, unsigned long previous);

// NEW: Function to adjust Kalman filter R matrix based on unbalance
void adjustKalmanR(bool unbalanced);
// NEW: Function to assess overall system health and potential unbalance
bool detectUnbalance(const PhaseAnalysisResult& result, float a_raw_offset, float b_raw_offset, float c_raw_offset);


unsigned long safeMicrosDiff(unsigned long current, unsigned long previous) {
    return (current >= previous) ? (current - previous) : ((ULONG_MAX - previous) + current + 1);
}

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println("Initializing Corrected BLDC Sensorless Encoder...");

    initM(); // Initialize Kalman matrices Qk1, Rk1, Ik1
    analogReference(INTERNAL);
    // Set ADC prescaler to 32 (16MHz / 32 = 500kHz). This is a good balance for 10-bit ADC.
    // For ATmega328P, ADPS2:ADPS0 bits are at ADCSRA[2:0].
    ADCSRA = (ADCSRA & ~0x07) | 0x05; // 0x05 corresponds to a prescaler of 32 (16MHz / 32 = 500kHz)
    // TCCR1B setup (for Timer1, not directly related to ADC, kept as in original)
    TCCR1B = (TCCR1B & 0b11111000) | 0x01; // No prescaling for Timer1 (CS12=0, CS11=0, CS10=1)

    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW); // Ensure CP is low initially
    aTEP();

    // Initialize ADC Kalman Filter state and covariance
    xk1.Fill(static_cast<float>(Config::TQA)); // Initialize state to expected quiescent ADC reading
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 100.0f; // High initial covariance for quick convergence

    // Initialize timestamps and angle for velocity calculation
    pVTs = micros();
    last_velocity_update_timestamp = micros();
    last_unwrapped_angle_for_velocity_calc = 0.0f;

    // Initial state for angle accumulation
    unwrappedMeasuredAngle = 0.0f; // Start at 0, or could be a known home position if available.
    lastRawElectricalAngle = 0.0f; // Align for initial updateUnwrappedMeasuredAngle calls.
    last_known_unwrapped_angle = 0.0f; // Initially no known history.
    estimatedVelocity = 0.0f; // Initialize estimated velocity

    Serial.println("System ready. Current state: UNLOCKED");
}

void loop() {
    unsigned long now = micros();
    unsigned long dt_micros = safeMicrosDiff(now, last_velocity_update_timestamp);
    float dt_sec = dt_micros / 1e6f;

    // Phase analysis includes DQ0 components and unbalance indicators
    PhaseAnalysisResult current_reading = analyzePhases(); 
    eSC(current_reading.magnitude); // Pass magnitude to gDNS for dynamic oversampling

    // Adjust Kalman filter's measurement noise based on detected unbalance
    adjustKalmanR(current_reading.is_unbalanced);

    // --- Update Unwrapped Angle & Calculate Instantaneous Velocity ---
    updateUnwrappedMeasuredAngle(current_reading.angle);

    float instantaneous_velocity = 0.0f;
    if (dt_sec > 1e-5) { // Prevent division by zero
        instantaneous_velocity = (unwrappedMeasuredAngle - last_unwrapped_angle_for_velocity_calc) / dt_sec;
    }

    // --- Smooth Estimated Velocity using EMA ---
    // Make smoothing more aggressive (lower alpha) if signal is unreliable or unbalanced
    float current_velocity_alpha = Config::VELOCITY_ALPHA;
    if (!current_reading.is_reliable || current_reading.is_unbalanced) {
        current_velocity_alpha *= 0.5f; // Halve alpha, means more smoothing / less responsiveness
    }

    if (current_reading.is_reliable) { // Only smooth if we have a reliable reading
        estimatedVelocity = current_velocity_alpha * instantaneous_velocity + (1.0f - current_velocity_alpha) * estimatedVelocity;
    } else {
        // If not reliable, slowly decay velocity to zero
        estimatedVelocity *= (1.0f - Config::CONFIDENCE_DECAY_RATE); // A small decay to zero
        if (fabsf(estimatedVelocity) < 0.001f) estimatedVelocity = 0.0f; // Snap to zero if very small
    }

    // --- Update `trackingConfidence` ---
    if (current_reading.is_reliable && fabsf(estimatedVelocity) > Config::VELOCITY_ENGAGE_THRESHOLD && !current_reading.is_unbalanced) {
        // Increase confidence if signal is reliable AND there's significant smoothed motion AND not unbalanced
        trackingConfidence += Config::CONFIDENCE_GAIN_RATE;
    } else if (!current_reading.is_reliable || fabsf(estimatedVelocity) < Config::VELOCITY_DISENGAGE_THRESHOLD || current_reading.is_unbalanced) {
        // Decrease confidence if signal is unreliable OR smoothed motion is negligible OR unbalanced
        trackingConfidence -= Config::CONFIDENCE_DECAY_RATE;
    }
    trackingConfidence = constrain(trackingConfidence, 0.0f, 1.0f); // Keep confidence between 0 and 1

    // --- Core State Machine (REVISED) ---
    switch (currentState) {
        case UNLOCKED: {
            estimatedAngle = unwrappedMeasuredAngle; // Keep estimated angle current with raw angle
            // Only attempt to engage if confidence is high AND not currently unbalanced
            if (trackingConfidence >= Config::ENGAGE_CONFIDENCE_THRESHOLD && !current_reading.is_unbalanced) {
                Serial.println("UNLOCKED -> ENGAGING: Confidence reached engage threshold.");
                // Re-synchronize unwrappedMeasuredAngle, accounting for missed revolutions
                if (has_been_locked_before) {
                    float currentWrappedAngle = current_reading.angle;
                    float lastKnownWrappedAngle = fmodf(last_known_unwrapped_angle, Config::TPF);
                    if (lastKnownWrappedAngle < 0) lastKnownWrappedAngle += Config::TPF;
                    if (currentWrappedAngle < 0) currentWrappedAngle += Config::TPF;

                    float angleDifference = currentWrappedAngle - lastKnownWrappedAngle;
                    if (angleDifference > Config::PF) angleDifference -= Config::TPF;
                    else if (angleDifference < -Config::PF) angleDifference += Config::TPF;

                    unwrappedMeasuredAngle = last_known_unwrapped_angle + angleDifference;
                } else {
                    unwrappedMeasuredAngle = current_reading.angle;
                }
                lastRawElectricalAngle = current_reading.angle; // Initialize for updateUnwrappedMeasuredAngle
                last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle; // Align velocity calc
                last_velocity_update_timestamp = now;

                currentState = ENGAGING;
            }
            break;
        }

        case ENGAGING: {
            if (trackingConfidence >= Config::LOCK_CONFIDENCE_THRESHOLD && !current_reading.is_unbalanced) {
                Serial.println("ENGAGING -> LOCKED: Confidence reached lock threshold.");
                currentState = LOCKED;
                has_been_locked_before = true; // We've successfully locked at least once.
            } else if (trackingConfidence < Config::UNLOCK_CONFIDENCE_THRESHOLD || current_reading.is_unbalanced) {
                // Drop to UNLOCKED if confidence is too low OR system becomes unbalanced
                Serial.print("ENGAGING -> UNLOCKED: ");
                if (trackingConfidence < Config::UNLOCK_CONFIDENCE_THRESHOLD) Serial.println("Confidence dropped below unlock threshold.");
                else Serial.println("Unbalance detected.");
                
                currentState = UNLOCKED;
                last_known_unwrapped_angle = unwrappedMeasuredAngle; // Store for re-sync
                estimatedVelocity = 0.0f; // Reset velocity
            }
            estimatedAngle = unwrappedMeasuredAngle; // Update estimated angle
            break;
        }

        case LOCKED: {
            if (trackingConfidence < Config::UNLOCK_CONFIDENCE_THRESHOLD || current_reading.is_unbalanced) {
                // Drop to UNLOCKED if confidence is too low OR system becomes unbalanced
                Serial.print("LOCKED -> UNLOCKED: ");
                if (trackingConfidence < Config::UNLOCK_CONFIDENCE_THRESHOLD) Serial.println("Confidence dropped below unlock threshold.");
                else Serial.println("Unbalance detected.");

                currentState = UNLOCKED;
                last_known_unwrapped_angle = unwrappedMeasuredAngle; // Store for re-sync
                estimatedVelocity = 0.0f; // Force velocity to zero upon losing lock.
            }
            estimatedAngle = unwrappedMeasuredAngle; // The most accurate angle when locked.
            break;
        }
    }

    // Update velocity tracking point for the next iteration. This must happen *after* state logic.
    last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle;
    last_velocity_update_timestamp = now;


    // Update final derived values for output
    eMS = round(estimatedAngle / Config::TPF * Config::ESR);
    rC = round(estimatedAngle / Config::TPF);

    pVTs = now;
    lD(current_reading); // Log data to Serial Plotter

    if (Serial.available() > 0) {
        char cmd = toupper(Serial.read());
        if (cmd == 'T') {
            aTEP(); // Call auto-tune if 'T' is received
        }
    }
}

// --- Utility and Support Functions ---

PhaseAnalysisResult analyzePhases() {
    PhaseAnalysisResult result;
    // Get the current Kalman filtered raw phase readings
    float a = xk1(0);
    float b = xk1(1);
    float c = xk1(2);

    // Subtract the common mode offset (Config::TQA).
    float a_offset_comp = a - Config::TQA;
    float b_offset_comp = b - Config::TQA;
    float c_offset_comp = c - Config::TQA;

    // --- Clarke Transform (Alpha-Beta components) ---
    float alpha = (2.0f / 3.0f) * (a_offset_comp - 0.5f * (b_offset_comp + c_offset_comp));
    float beta  = (2.0f / 3.0f) * (Config::SQ3 / 2.0f) * (b_offset_comp - c_offset_comp);

    // Calculate angle from Clarke components
    result.angle = atan2f(beta, alpha); // angle will be in (-PI, PI]
    result.magnitude = sqrtf(alpha * alpha + beta * beta);
    result.is_reliable = result.magnitude > Config::MIN_SIGNAL_MAGNITUDE;

    // --- Park Transform (DQ components) ---
    float cos_theta_e = cosf(result.angle);
    float sin_theta_e = sinf(result.angle);

    result.d_component = alpha * cos_theta_e + beta * sin_theta_e;
    result.q_component = -alpha * sin_theta_e + beta * cos_theta_e;

    // --- Zero-sequence component (0 component) ---
    result.zero_component = (a_offset_comp + b_offset_comp + c_offset_comp) / 3.0f;

    // --- Unbalance Detection --- (NEW LOGIC)
    result.is_unbalanced = detectUnbalance(result, a_offset_comp, b_offset_comp, c_offset_comp);

    // Store Kalman filter residual sum for logging/analysis
    // Note: The actual residual calculation happens in aMKF, this is just for logging the current iteration's residual.
    // To get the sum of absolute residuals here, you'd need access to the `y` matrix from the last `aMKF` call.
    // For simplicity, let's assume we can get it or re-calculate it for logging.
    // A better approach would be to pass `y` out of `aMKF` or make it a global for analysis.
    // For now, let's make the residual sum calculation part of `aMKF` and store it in a global.
    // This is a placeholder for where that sum would be used to set `result.kalman_residual_sum_abs`.
    // It's more accurate to report the actual residual from the Kalman update step.
    // For now, let's keep it simple and approximate, or assume it's set by `aMKF`.

    return result;
}

// NEW FUNCTION: Detect unbalance based on various metrics
bool detectUnbalance(const PhaseAnalysisResult& result, float a_offset_comp, float b_offset_comp, float c_offset_comp) {
    // 1. Check Zero-sequence component
    if (fabsf(result.zero_component) > Config::ZERO_COMPONENT_UNBALANCE_THRESHOLD) {
        // Serial.println("Unbalance detected: High zero component"); // For debugging
        return true;
    }

    // 2. Check magnitude balance of individual phases (after offset compensation)
    // Only meaningful if there's significant signal magnitude overall
    if (result.magnitude > Config::MIN_SIGNAL_MAGNITUDE * 2) { // Only check if magnitude is reasonable
        float max_phase_val = fmaxf(fabsf(a_offset_comp), fmaxf(fabsf(b_offset_comp), fabsf(c_offset_comp)));
        float min_phase_val = fminf(fabsf(a_offset_comp), fminf(fabsf(b_offset_comp), fabsf(c_offset_comp)));
        
        if (max_phase_val > 0 && (max_phase_val - min_phase_val) / max_phase_val > Config::PHASE_MAGNITUDE_IMBALANCE_RATIO) {
            // Serial.println("Unbalance detected: Phase magnitude imbalance"); // For debugging
            return true;
        }
    }
    
    // 3. Check D/Q component behavior, especially during low velocity/stationary periods
    // If estimated velocity is very low, D and Q should be relatively stable/small (if BEMF-based sensing)
    // This assumes the D-axis is aligned with the flux and Q with torque.
    if (fabsf(estimatedVelocity) < Config::VELOCITY_ENGAGE_THRESHOLD) { // Motor is nearly stationary or just starting
        if (fabsf(result.d_component) > Config::DQ_UNBALANCE_THRESHOLD || fabsf(result.q_component) > Config::DQ_UNBALANCE_THRESHOLD) {
            // This threshold needs tuning. If d_component should ideally be 0, a non-zero value indicates misalignment/unbalance.
            // If q_component is high but velocity is low, it suggests an issue.
            // Serial.println("Unbalance detected: High D/Q at low speed"); // For debugging
            return true;
        }
    }

    // 4. Check Kalman Filter residuals (from previous update, as it's updated *before* this call in loop)
    // We need a way to get the *previous* residual sum from aMKF.
    // For now, let's make `kalman_residual_sum_abs` a global or pass it from `aMKF`.
    // For this example, let's modify `aMKF` to update a global variable `last_kalman_residual_sum_abs`.
    // Then `analyzePhases` can pick it up.

    if (result.kalman_residual_sum_abs > Config::KALMAN_RESIDUAL_THRESHOLD) {
         // Serial.println("Unbalance detected: High Kalman residual"); // For debugging
        return true;
    }


    return false; // No significant unbalance detected
}

void updateUnwrappedMeasuredAngle(float currentRawAngle) {
    float angleDelta = currentRawAngle - lastRawElectricalAngle;
    // Normalize angleDelta to be within (-PI, PI)
    if (angleDelta > Config::PF) angleDelta -= Config::TPF;
    else if (angleDelta < -Config::PF) angleDelta += Config::TPF;

    unwrappedMeasuredAngle += angleDelta;
    lastRawElectricalAngle = currentRawAngle;
}

void initM() {
    Ik1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Ik1(i, i) = 1.0f;
    Qk1.Fill(0.0f); 
    Rk1.Fill(0.0f); // Will be set by adjustKalmanR initially and dynamically
    // Initialize Rk1 with base value
    for (int i = 0; i < NKF1; i++) {
        Qk1(i, i) = Config::QDK1_BASE;
        Rk1(i, i) = Config::RDK1_BASE; 
    }
}

// NEW FUNCTION: Adjust Kalman filter R matrix
void adjustKalmanR(bool unbalanced) {
    // Dynamically adjust R (measurement noise covariance) based on unbalance
    if (unbalanced) {
        // Increase R to trust measurements less, rely more on model
        for (int i = 0; i < NKF1; i++) {
            Rk1(i, i) = Config::RDK1_BASE * 5.0f; // Example: 5x higher measurement noise
        }
    } else {
        // Revert to base R when balanced
        for (int i = 0; i < NKF1; i++) {
            Rk1(i, i) = Config::RDK1_BASE;
        }
    }
    // You could also add dynamic adjustment for Q (process noise) if motor dynamics are explicitly modeled and change.
    // For sensorless BEMF, usually R is varied more.
}


void aTEP() {
    Serial.println("\n--- Starting Auto-Tune (Ensure motor is stationary) ---");
    xk1.Fill(0.0f); Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 1000.0f;
    for (int iter = 0; iter < Config::PCI; iter++) {
        long rT[3] = {0};
        for (int i = 0; i < Config::DNS; i++) {
            // Original excitation pulse
            pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW); TCNT1 = 0;
            analogWrite(Config::CP, EPW); delayMicroseconds(Config::EPWU);
            analogWrite(Config::CP, 0); // Turn off excitation

            // Dithering occurs within cAR now for each individual ADC read
            rT[0] += cAR(Config::PA);
            rT[1] += cAR(Config::PB);
            rT[2] += cAR(Config::PC);
        }
        BLA::Matrix<NKF1, 1> avg_rT;
        for (int i = 0; i < 3; i++) avg_rT(i) = static_cast<float>(rT[i]) / static_cast<float>(Config::DNS);
        aMKF(avg_rT); // Update Kalman filter
        int cAA = static_cast<int>((xk1(0) + xk1(1) + xk1(2)) / 3.0f);
        int err = cAA - Config::TQA;
        if (abs(err) <= Config::PCT) { Serial.println("--- Auto-tune complete! ---"); break; }
        EPW += (err < 0) ? Config::PAS : -Config::PAS;
        EPW = constrain(EPW, 1, 255);
        //delay(5); // Consider removing or using minimal delay for faster tuning
        Serial.print("err:"); Serial.print(err);
        Serial.print(",PWM :"); Serial.println(EPW);
    }
    Serial.print("Final Tuned PWM Value: "); Serial.println(EPW);
}

// Modified: gDNS now takes signal magnitude
int gDNS(float currentMagnitude) {
    if (currentMagnitude < Config::LOW_MAGNITUDE_THRESHOLD) {
        return Config::MAX_S; // More samples for weak signals
    }
    if (currentMagnitude > Config::HIGH_MAGNITUDE_THRESHOLD) {
        return Config::MIN_S; // Fewer samples for strong signals
    }
    // Linearly interpolate between MAX_S and MIN_S based on magnitude
    float normalizedMagnitude = (currentMagnitude - Config::LOW_MAGNITUDE_THRESHOLD) /
                                (Config::HIGH_MAGNITUDE_THRESHOLD - Config::LOW_MAGNITUDE_THRESHOLD);
    int numSamples = Config::MAX_S - static_cast<int>(normalizedMagnitude * (Config::MAX_S - Config::MIN_S));
    return constrain(numSamples, Config::MIN_S, Config::MAX_S);
}

// Modified: eSC now takes current_signal_magnitude
void eSC(float current_signal_magnitude) {
    long r[3] = {0};
    int nST = gDNS(current_signal_magnitude); // Pass magnitude to gDNS
        pinMode(Config::CP, OUTPUT);
        digitalWrite(Config::CP, LOW); // Ensure CP is low for a clean start
        TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.
        analogWrite(Config::CP, EPW);
        delayMicroseconds(Config::EPWU);

    for (int i = 0; i < nST; i++) {
        r[0] += cAR(Config::PA);
        r[1] += cAR(Config::PB);
        r[2] += cAR(Config::PC);
    }
    BLA::Matrix<NKF1, 1> avg_r;
    for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(r[i]) / static_cast<float>(nST);
    aMKF(avg_r); // Update Kalman filter with new measurements
    TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.
    analogWrite(Config::CP, EPW); // Keep excitation pulse active after readings
}

/**
 * @brief Reads an analog pin with dithering applied.
 * This function applies a short pulse on the CP pin before taking an ADC reading
 * to introduce dithering noise, improving ADC linearity and effective resolution
 * when combined with oversampling.
 * @param p The analog pin to read (e.g., A0, A1, A2).
 * @return The analog reading.
 */
int cAR(int p) {
    // Apply a dither pulse before conversion
    TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.
    analogWrite(Config::CP, Config::DITHER_AMPLITUDE_PIN_STATE);
    delayMicroseconds(Config::DITHER_PULSE_WIDTH_US);
    TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.

    analogWrite(Config::CP, EPW); // Restore main excitation PWM
    return analogRead(p);
}

// Global variable to store the last Kalman residual sum for analysis
float last_kalman_residual_sum_abs_global = 0.0f;

void aMKF(const BLA::Matrix<NKF1, 1>& z) {
    Pk1 += Qk1; // Prediction step (covariance)

    BLA::Matrix<NKF1, 1> y = z - xk1; // Innovation (measurement residual)

    // Calculate sum of absolute residuals for unbalance detection
    float current_residual_sum_abs = 0.0f;
    for(int i = 0; i < NKF1; ++i) {
        current_residual_sum_abs += fabsf(y(i));
    }
    last_kalman_residual_sum_abs_global = current_residual_sum_abs; // Store for `analyzePhases`

    BLA::Matrix<NKF1, NKF1> S = Pk1 + Rk1; // Innovation covariance
    
    // Check for singularity before inverting
    float det_S = BLA::Determinant(S);
    if (fabsf(det_S) < Config::MATRIX_SINGULARITY_THRESHOLD) {
        // Serial.println("Kalman: S matrix close to singular, skipping update.");
        return; // Skip update if S is singular to prevent errors
    }
    auto K = Pk1 * BLA::Inverse(S); // Kalman Gain

    xk1 += K * y; // Update step (state estimate)
    Pk1 = (Ik1 - K) * Pk1; // Update step (covariance)
}

void lD(const PhaseAnalysisResult& result) {
    // For Arduino IDE Serial Plotter
    Serial.print("A:"); Serial.print(xk1(0)-Config::TQA); Serial.print(",");
    Serial.print("B:"); Serial.print(xk1(1)-Config::TQA); Serial.print(",");
    Serial.print("C:"); Serial.print(xk1(2)-Config::TQA); Serial.print(",");
    Serial.print("Angle:"); Serial.print(estimatedAngle, 2); Serial.print(",");
    Serial.print("Vel:"); Serial.print(estimatedVelocity, 3); Serial.print(",");
    Serial.print("Mag:"); Serial.print(result.magnitude, 2); Serial.print(",");
    Serial.print("Conf:"); Serial.print(trackingConfidence, 2); Serial.print(",");
    Serial.print("D:"); Serial.print(result.d_component, 2); Serial.print(",");
    Serial.print("Q:"); Serial.print(result.q_component, 2); Serial.print(",");
    Serial.print("0:"); Serial.print(result.zero_component, 2); Serial.print(",");
    Serial.print("Unbal:"); Serial.print(result.is_unbalanced ? 1 : 0); Serial.print(","); // Log unbalance flag
    Serial.print("KRes:"); Serial.print(last_kalman_residual_sum_abs_global, 2); Serial.print(","); // Log Kalman residual sum
    Serial.print("eMS:"); Serial.print(eMS);Serial.print(",");
    Serial.print("rC:"); Serial.print(rC); Serial.print(",");
    Serial.print("State:"); Serial.println(currentState); // 0=UNLOCKED, 1=ENGAGING, 2=LOCKED
}
