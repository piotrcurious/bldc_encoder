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
    constexpr int MAX_S = 128; // MAX_NUM_SAMPLES

    // NEW: Magnitude Thresholds for Dynamic Oversampling
    constexpr float LOW_MAGNITUDE_THRESHOLD = 5.0f; // Below this, use MAX_S samples
    constexpr float HIGH_MAGNITUDE_THRESHOLD = 7.0f; // Above this, use MIN_S samples
                                                    // Adjust these values based on your ADC's typical magnitude range

    // Dithering Configuration
    constexpr int DITHER_PULSE_WIDTH_US = 32; // Duration of the dither pulse in microseconds
                                            // This might need fine-tuning.
    constexpr int DITHER_AMPLITUDE_PIN_STATE = 32; // Or LOW, or an analogWrite value for more subtle dither

    // Excitation & Auto-Tuning
    constexpr int EPWU = 256, FSTU = 1;
    constexpr int TQA = 32, PAS = 1, PCI = 2000, PCT = 1;

    // Kalman Filter (for ADC smoothing)
    constexpr float QDK1 = 20.0f, RDK1 = 20.1f;
    constexpr float MATRIX_SINGULARITY_THRESHOLD = 1e-3f;

    // --- Core Estimation & Tracking Parameters ---
    constexpr float MIN_SIGNAL_MAGNITUDE = 0.05f; // Minimum vector magnitude for a reliable signal (in ADC units). This is a critical tuning parameter.
    constexpr float VST = 0.01f; // Minimum velocity to consider moving (rad/s).
    constexpr int ENGAGE_CONFIRMATION_FRAMES = 16; // How many frames of consistent motion are needed to lock tracking.
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
int consistentMovementCounter = 0;

// Angle & Velocity State
float lastRawElectricalAngle = 0.0f; // Stores the raw angle from the *previous* reading, used by updateUnwrappedMeasuredAngle.
float unwrappedMeasuredAngle = 0.0f; // The continuously tracked electrical angle.
float estimatedAngle = 0.0f; // Output angle, usually unwrappedMeasuredAngle when locked.
float estimatedVelocity = 0.0f; // Output velocity.

// NEW: Store the unwrapped angle when we lose lock
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

unsigned long safeMicrosDiff(unsigned long current, unsigned long previous) {
    return (current >= previous) ? (current - previous) : ((ULONG_MAX - previous) + current + 1);
}

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println("Initializing Corrected BLDC Sensorless Encoder...");

    initM();
    analogReference(INTERNAL);
    // Set ADC prescaler to 32 (16MHz / 32 = 500kHz). This is a good balance for 10-bit ADC.
    // For ATmega328P, ADPS2:ADPS0 bits are at ADCSRA[2:0].
    ADCSRA = (ADCSRA & ~0x07) | 0x05; // 0x05 corresponds to a prescaler of 32 (16MHz / 32 = 500kHz)
    // TCCR1B setup (for Timer1, not directly related to ADC, kept as in original)
    TCCR1B = (TCCR1B & 0b11111000) | 0x01; // No prescaling for Timer1 (CS12=0, CS11=0, CS10=1)

    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW); // Ensure CP is low initially
    aTEP();

    // Initialize ADC Kalman Filter
    xk1.Fill(static_cast<float>(Config::TQA));
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 100.0f;

    // Initialize timestamps and angle for velocity calculation
    pVTs = micros();
    last_velocity_update_timestamp = micros();
    last_unwrapped_angle_for_velocity_calc = 0.0f;

    // Initial state for angle accumulation
    unwrappedMeasuredAngle = 0.0f; // Start at 0, or could be a known home position if available.
    lastRawElectricalAngle = 0.0f; // Align for initial updateUnwrappedMeasuredAngle calls.
    last_known_unwrapped_angle = 0.0f; // Initially no known history.

    Serial.println("System ready. Current state: UNLOCKED");
}

void loop() {
    unsigned long now = micros();
    unsigned long dt_micros = safeMicrosDiff(now, last_velocity_update_timestamp);
    float dt_sec = dt_micros / 1e6f;

    PhaseAnalysisResult current_reading = analyzePhases(); // Get angle and signal quality
    eSC(current_reading.magnitude); // Modified: Pass magnitude to eSC

    // --- Core State Machine ---
    switch (currentState) {
        case UNLOCKED: {
            estimatedVelocity = 0.0f;
            consistentMovementCounter = 0; // Reset counter in UNLOCKED state

            if (current_reading.is_reliable) {
                Serial.println("UNLOCKED -> ENGAGING: Reliable signal detected.");

                // CRITICAL CHANGE: Re-synchronize unwrappedMeasuredAngle without resetting accumulated value
                if (has_been_locked_before) {
                    // Calculate the difference between the current raw angle and the wrapped
                    // portion of the last known unwrapped angle.
                    float currentWrappedAngle = current_reading.angle;
                    float lastKnownWrappedAngle = fmodf(last_known_unwrapped_angle, Config::TPF);
                    if (lastKnownWrappedAngle < 0) { // Ensure fmodf result is always [0, TPF) or (-TPF, 0] depending on fmodf
                        lastKnownWrappedAngle += Config::TPF; // Normalize to [0, TPF)
                    }
                    if (currentWrappedAngle < 0) {
                        currentWrappedAngle += Config::TPF; // Normalize to [0, TPF)
                    }

                    // Calculate the shortest angular distance between current and last known wrapped angle.
                    float angleDifference = currentWrappedAngle - lastKnownWrappedAngle;
                    if (angleDifference > Config::PF) {
                        angleDifference -= Config::TPF;
                    } else if (angleDifference < -Config::PF) {
                        angleDifference += Config::TPF;
                    }

                    // Add this difference to the *last known unwrapped angle* to get the new starting point
                    // for unwrappedMeasuredAngle. This accounts for full revolutions missed.
                    unwrappedMeasuredAngle = last_known_unwrapped_angle + angleDifference;

                } else {
                    // First time engaging, start unwrapped angle from current raw reading.
                    unwrappedMeasuredAngle = current_reading.angle;
                }

                lastRawElectricalAngle = current_reading.angle; // Initialize for updateUnwrappedMeasuredAngle

                // Reset the velocity tracking point to the current unwrapped angle
                // and current time for this new engagement.
                last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle;
                last_velocity_update_timestamp = now;

                currentState = ENGAGING;
            }
            break;
        }

        case ENGAGING: {
            estimatedVelocity = 0.0f; // Still output zero velocity while engaging

            if (!current_reading.is_reliable) {
                Serial.println("ENGAGING -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                consistentMovementCounter = 0;
                // When losing signal, store the current unwrapped angle
                // for potential re-synchronization later.
                last_known_unwrapped_angle = unwrappedMeasuredAngle;
                break;
            }

            // Update unwrappedMeasuredAngle based on current and last raw readings.
            // This allows `unwrappedMeasuredAngle` to continue tracking even in ENGAGING.
            updateUnwrappedMeasuredAngle(current_reading.angle);

            float instant_vel = 0.0f;
            if (dt_sec > 1e-5) {
                instant_vel = (unwrappedMeasuredAngle - last_unwrapped_angle_for_velocity_calc) / dt_sec;
            }

            // Update velocity tracking point for the next iteration in ENGAGING.
            last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle;
            last_velocity_update_timestamp = now;

            // Check for consistent, non-zero movement.
            if (fabsf(instant_vel) > Config::VST) {
                consistentMovementCounter++;
            } else {
                consistentMovementCounter = 0; // Reset if movement stops or is negligible.
            }

            if (consistentMovementCounter >= Config::ENGAGE_CONFIRMATION_FRAMES && current_reading.is_reliable) {
                Serial.println("ENGAGING -> LOCKED: Motion confirmed.");
                currentState = LOCKED;
                has_been_locked_before = true; // We've successfully locked at least once.

                // Align velocity tracking for LOCKED state.
                last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle;
                last_velocity_update_timestamp = now;

                estimatedAngle = unwrappedMeasuredAngle;
                estimatedVelocity = instant_vel; // Seed with the last instantaneous velocity.
            }
            break;
        }

        case LOCKED: {
            if (!current_reading.is_reliable) {
                Serial.println("LOCKED -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                estimatedVelocity = 0.0f; // Force velocity to zero.
                // Store the current unwrapped angle for potential re-synchronization later.
                last_known_unwrapped_angle = unwrappedMeasuredAngle;
                break;
            }

            // First, update the unwrapped angle for the current reading
            updateUnwrappedMeasuredAngle(current_reading.angle);

            // Now, calculate velocity using the newly updated unwrappedMeasuredAngle
            // and the values from the *previous* loop iteration.
            float instantaneous_velocity = 0.0f;
            if (dt_sec > 1e-5) {
                instantaneous_velocity = (unwrappedMeasuredAngle - last_unwrapped_angle_for_velocity_calc) / dt_sec;
            } else {
                instantaneous_velocity = 0.0f;
            }

            // Update `last_unwrapped_angle_for_velocity_calc` and `last_velocity_update_timestamp`
            // for the *next* loop iteration's velocity calculation.
            last_unwrapped_angle_for_velocity_calc = unwrappedMeasuredAngle;
            last_velocity_update_timestamp = now;

            estimatedVelocity = instantaneous_velocity;
            estimatedAngle = unwrappedMeasuredAngle; // The most accurate angle when locked.
            break;
        }
    }

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
    float a = xk1(0), b = xk1(1), c = xk1(2);

    a -= Config::TQA; b -= Config::TQA; c -= Config::TQA;

    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta  = (2.0f / 3.0f) * (Config::SQ3 / 2.0f) * (b - c);

    result.angle = atan2f(beta, alpha); // angle will be in (-PI, PI]
    result.magnitude = sqrtf(alpha * alpha + beta * beta);
    result.is_reliable = result.magnitude > Config::MIN_SIGNAL_MAGNITUDE;
    return result;
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
    Qk1.Fill(0.0f); Rk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) {
        Qk1(i, i) = Config::QDK1;
        Rk1(i, i) = Config::RDK1;
    }
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
        aMKF(avg_rT);
        int cAA = static_cast<int>((xk1(0) + xk1(1) + xk1(2)) / 3.0f);
        int err = cAA - Config::TQA;
        if (abs(err) <= Config::PCT) { Serial.println("--- Auto-tune complete! ---"); break; }
        EPW += (err < 0) ? Config::PAS : -Config::PAS;
        EPW = constrain(EPW, 1, 255);
        //delay(5);
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
        // Original excitation pulse
        // Note: The original code sets CP to LOW and then analogWrites EPW.
        // For sensing, the CP pin should typically be high-impedance (INPUT)
        // or driven low *after* the excitation pulse, and then possibly dithered.
        // I'm keeping the original excitation pattern here but ensuring it finishes
        // before the individual ADC reads within cAR.
//        pinMode(Config::CP, OUTPUT);
//        digitalWrite(Config::CP, LOW); // Ensure CP is low for a clean start
//        TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.
//        analogWrite(Config::CP, EPW);
//        delayMicroseconds(Config::EPWU);
        //analogWrite(Config::CP, 0); // Turn off the excitation after pulse width

        // Dithering happens inside cAR for each individual reading
        r[0] += cAR(Config::PA);
        r[1] += cAR(Config::PB);
        r[2] += cAR(Config::PC);
    }
    BLA::Matrix<NKF1, 1> avg_r;
    for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(r[i]) / static_cast<float>(nST);
    aMKF(avg_r);

    // After all samples are taken for the current measurement cycle,
    // you might want to return CP to its idle state (e.g., HIGH for charge pump)
    // or keep it low if that's its typical resting state outside of excitation.
    // The original code had `digitalWrite(Config::CP, EPW);` here which is unusual
    // since EPW is an analogWrite value. Assuming it means `digitalWrite(Config::CP, HIGH);`
    // or similar if CP needs to stay active for some reason. For now, setting it to low.
    //digitalWrite(Config::CP, LOW);
    TCNT1 = 0; // Reset Timer1, presumably for timing the excitation pulse accurately.
    analogWrite(Config::CP, EPW);
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

    //digitalWrite(Config::CP, LOW); // Return CP to a stable state (e.g., LOW or its default idle)
    analogWrite(Config::CP, EPW);
    return analogRead(p);
}


void aMKF(const BLA::Matrix<NKF1, 1>& z) {
    Pk1 += Qk1;
    BLA::Matrix<NKF1, 1> y = z - xk1;
    BLA::Matrix<NKF1, NKF1> S = Pk1 + Rk1;
    if (fabsf(BLA::Determinant(S)) < Config::MATRIX_SINGULARITY_THRESHOLD) return;
    auto K = Pk1 * BLA::Inverse(S);
    xk1 += K * y;
    Pk1 = (Ik1 - K) * Pk1;
}


void lD(const PhaseAnalysisResult& result) {
    // For Arduino IDE Serial Plotter
    Serial.print(xk1(0)-Config::TQA); Serial.print(",");
    Serial.print(xk1(1)-Config::TQA); Serial.print(",");
    Serial.print(xk1(2)-Config::TQA); Serial.print(",");
    Serial.print("Angle:"); Serial.print(estimatedAngle, 2); Serial.print(",");
    Serial.print("Vel:"); Serial.print(estimatedVelocity, 3); Serial.print(",");
    Serial.print("Mag:"); Serial.print(result.magnitude, 2); Serial.print(",");
    Serial.print(eMS);Serial.print(",");
    Serial.print(rC);
    Serial.print("State:"); Serial.println(currentState); // 0=UNLOCKED, 1=ENGAGING, 2=LOCKED
}
