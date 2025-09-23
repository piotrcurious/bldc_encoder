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
    constexpr int DNS = 8; // DEFAULT_NUM_SAMPLES
    constexpr int MIN_S = 4; // MIN_NUM_SAMPLES
    constexpr int MAX_S = 8; // MAX_NUM_SAMPLES
    constexpr float LV_TH = 5.0f; // LOW_VELOCITY_THRESHOLD (rad/s)
    constexpr float HV_TH = 15.0f; // HIGH_VELOCITY_THRESHOLD (rad/s)

    // Excitation & Auto-Tuning
    constexpr int EPWU = 512, FSTU = 1;
    constexpr int TQA = 32, PAS = 1, PCI = 2000, PCT = 1;

    // Kalman Filter (for ADC smoothing)
    constexpr float QDK1 = 1.0f, RDK1 = 80.1f;
    constexpr float MATRIX_SINGULARITY_THRESHOLD = 1e-5f;

    // --- Core Estimation & Tracking Parameters ---
    constexpr float MIN_SIGNAL_MAGNITUDE = 1.0f; // Minimum vector magnitude for a reliable signal (in ADC units). This is a critical tuning parameter.
    constexpr float VST = 0.01f; // Minimum velocity to consider moving (rad/s).
    constexpr int ENGAGE_CONFIRMATION_FRAMES = 10; // How many frames of consistent motion are needed to lock tracking.
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
float lastRawElectricalAngle = 0.0f;
float unwrappedMeasuredAngle = 0.0f;
float estimatedAngle = 0.0f;
float estimatedVelocity = 0.0f;

// Timestamps and values for velocity calculation
unsigned long pVTs = 0;
unsigned long last_update_timestamp = 0;
float last_unwrapped_angle_for_vel = 0.0f;

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
void eSC();
PhaseAnalysisResult analyzePhases();
void updateUnwrappedMeasuredAngle(float currentRawAngle);
void lD(const PhaseAnalysisResult& result);
int cAR(int);
void aMKF(const BLA::Matrix<NKF1, 1>& z);
int gDNS(float velocity);
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
    ADCSRA = (ADCSRA & ~0x07) | 0x05;
    TCCR1B = (TCCR1B & 0b11111000) | 0x01;

    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW);
    aTEP();

    // Initialize ADC Kalman Filter
    xk1.Fill(static_cast<float>(Config::TQA));
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 100.0f;
    
    // Initialize timestamps
    pVTs = micros();
    last_update_timestamp = pVTs;

    Serial.println("System ready. Current state: UNLOCKED");
}

void loop() {
    unsigned long now = micros();

    eSC(); // Perform ADC sampling cycle
    PhaseAnalysisResult reading = analyzePhases(); // Get angle and signal quality

    // --- Core State Machine ---
    switch (currentState) {
        case UNLOCKED: {
            estimatedVelocity = 0.0f;
            consistentMovementCounter = 0; // Reset counter in UNLOCKED state

            if (reading.is_reliable) {
                Serial.println("UNLOCKED -> ENGAGING: Reliable signal detected.");
                
                // --- CRITICAL CHANGE HERE ---
                // Re-synchronize unwrappedMeasuredAngle with the current reading.angle
                // to account for any un-tracked movement while UNLOCKED.
                float currentWrappedAngle = reading.angle; // The raw angle from -PI to PI
                float lastUnwrappedRemainder = fmodf(unwrappedMeasuredAngle, Config::TPF);
                if (lastUnwrappedRemainder < 0) { // Ensure fmodf result is always [0, TPF)
                    lastUnwrappedRemainder += Config::TPF;
                }
                
                // Calculate the difference between the current wrapped angle and
                // the wrapped part of the last unwrapped angle.
                float angleDiffWrapped = currentWrappedAngle - lastUnwrappedRemainder;

                // Adjust angleDiffWrapped to be in (-PI, PI)
                if (angleDiffWrapped > Config::PF) {
                    angleDiffWrapped -= Config::TPF;
                } else if (angleDiffWrapped < -Config::PF) {
                    angleDiffWrapped += Config::TPF;
                }

                // Add the adjusted wrapped difference to the current unwrapped angle
                // to effectively shift it by the un-tracked revolution(s).
                unwrappedMeasuredAngle += angleDiffWrapped;

                // Now, lastRawElectricalAngle can be set to the current reading.angle,
                // as unwrappedMeasuredAngle is correctly aligned.
                lastRawElectricalAngle = reading.angle;
                // --- END CRITICAL CHANGE ---
                
                last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
                last_update_timestamp = now;

                currentState = ENGAGING;
            }
            break;
        }
        
        case ENGAGING: {
            estimatedVelocity = 0.0f; // Output zero velocity while engaging

            if (!reading.is_reliable) {
                Serial.println("ENGAGING -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                consistentMovementCounter = 0;
                // DO NOT RESET unwrappedMeasuredAngle here. It should hold its last value.
                break;
            }

            // In ENGAGING state, we update `lastRawElectricalAngle` but *not* `unwrappedMeasuredAngle`
            // with `angleDelta` from `updateUnwrappedMeasuredAngle`.
            // Instead, we just keep `lastRawElectricalAngle` aligned with the current raw reading.
            // The `unwrappedMeasuredAngle` should only truly accumulate in the LOCKED state.
            // The critical part here is *not* calling `updateUnwrappedMeasuredAngle` which
            // directly modifies `unwrappedMeasuredAngle`.
            lastRawElectricalAngle = reading.angle; // Keep lastRawElectricalAngle updated for the next calculation.

            // Calculate instantaneous velocity to check for stable movement.
            // This calculation still needs to happen to determine if we can transition to LOCKED.
            // The velocity here is based on the *difference* between the current wrapped angle
            // and the previous wrapped angle, but it's not meant to contribute to `unwrappedMeasuredAngle` directly yet.
            float instant_vel = 0.0f;
            float dt_sec = safeMicrosDiff(now, last_update_timestamp) / 1e6f;
            if (dt_sec > 1e-5) {
                // Calculate instantaneous velocity based on the *change* in the current raw angle
                // relative to the last time we updated `last_unwrapped_angle_for_vel` (which was when entering ENGAGING).
                // This velocity is used for the state transition logic, not for accumulating angle.
                // You might need to refine how `instant_vel` is calculated in ENGAGING if `unwrappedMeasuredAngle`
                // isn't meant to be updated. For simplicity, let's keep it as is for the state transition,
                // as `unwrappedMeasuredAngle` was just resynchronized on entry to ENGAGING.
                instant_vel = (reading.angle - fmodf(last_unwrapped_angle_for_vel, Config::TPF)); // Approximate velocity for state transition check
                // Normalize to (-PI, PI)
                if (instant_vel > Config::PF) instant_vel -= Config::TPF;
                else if (instant_vel < -Config::PF) instant_vel += Config::TPF;
                instant_vel /= dt_sec;
            }
            
            last_update_timestamp = now;
            //last_unwrapped_angle_for_vel = unwrappedMeasuredAngle; // This line should *not* be here in ENGAGING if we want to avoid jumps

            // Check for consistent, non-zero movement.
            if (fabsf(instant_vel) > Config::VST) {
                consistentMovementCounter++;
            } else {
                consistentMovementCounter = 0; // Reset if movement stops or is negligible.
            }

            // Check if unwrapped angle matches last estimated unwrapped angle
            // This check might be redundant or less important with the re-synchronization above,
            // but keeping it for now if you found it useful.
            constexpr float ANGLE_ALIGNMENT_TOLERANCE = 0.01f; // Adjust as needed (e.g., 0.01 radians)
            // Note: The original anglesAligned calculation was a bit convoluted due to `fmodf`
            // potentially returning negative results. A simpler check might be:
            // fabsf(angleDelta) < ANGLE_ALIGNMENT_TOLERANCE
            // where angleDelta is the unwrapped delta from current `reading.angle` to `lastRawElectricalAngle`.
            // Given the updateUnwrappedMeasuredAngle logic, `lastRawElectricalAngle` is just the previous `reading.angle`.
            // The purpose of this `anglesAligned` check is to see if the current raw measurement
            // is "close enough" to where the tracking *expects* it to be.
            bool anglesAligned = fabsf(fmodf(unwrappedMeasuredAngle - reading.angle, Config::TPF)) < ANGLE_ALIGNMENT_TOLERANCE ||
                                 fabsf(fmodf(unwrappedMeasuredAngle - reading.angle, Config::TPF) - Config::TPF) < ANGLE_ALIGNMENT_TOLERANCE ||
                                 fabsf(fmodf(unwrappedMeasuredAngle - reading.angle, Config::TPF) + Config::TPF) < ANGLE_ALIGNMENT_TOLERANCE;
            // A more robust way to check alignment would be based on the instantaneous velocity
            // being consistent with the expected velocity, and the raw angle not jumping.
            // The `updateUnwrappedMeasuredAngle` already handles the wrapped nature for *small* steps.
            // The main purpose of this `anglesAligned` is to confirm the initial phase alignment.

            if (consistentMovementCounter > Config::ENGAGE_CONFIRMATION_FRAMES && reading.is_reliable /* && anglesAligned */) { // Removed anglesAligned from condition for initial testing
                Serial.println("ENGAGING -> LOCKED: Motion confirmed.");
                currentState = LOCKED;
                // When transitioning to LOCKED, initialize the `last_unwrapped_angle_for_vel`
                // and `estimatedAngle` to the current `unwrappedMeasuredAngle`
                // which was resynchronized when entering ENGAGING.
                last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
                estimatedAngle = unwrappedMeasuredAngle; // Initialize estimatedAngle with the currently unwrapped angle
                estimatedVelocity = instant_vel; // Give the filter a starting seed value.
            }
            break;
        }

        case LOCKED: {
            if (!reading.is_reliable) {
                Serial.println("LOCKED -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                estimatedVelocity = 0.0f; // Force velocity to zero.
                // consistentMovementCounter is reset when entering UNLOCKED.
                // DO NOT RESET unwrappedMeasuredAngle here. It must hold its value.
                break;
            }

            // ONLY update unwrappedMeasuredAngle with the change during the LOCKED state
            updateUnwrappedMeasuredAngle(reading.angle);
            
            // Calculate instantaneous velocity
            float instantaneous_velocity = 0.0f;
            float dt_sec = safeMicrosDiff(now, last_update_timestamp) / 1e6f;
            if (dt_sec > 1e-5) {
                instantaneous_velocity = (unwrappedMeasuredAngle - last_unwrapped_angle_for_vel) / dt_sec;
            }

            last_update_timestamp = now;
            last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
            
            // Directly assign instantaneous_velocity to estimatedVelocity
            estimatedVelocity = instantaneous_velocity;
            
            // In the LOCKED state, the estimated angle follows the reliable measurement.
            estimatedAngle = unwrappedMeasuredAngle;
            break;
        }
    }

    // Update final derived values for output
    eMS = round(estimatedAngle / Config::TPF * Config::ESR);
    rC = round(estimatedAngle / Config::TPF);
    
    pVTs = now;
    lD(reading); // Log data to Serial Plotter

    if (Serial.available() > 0) {
        char cmd = toupper(Serial.read());
        if (cmd == 'T') {
            aTEP(); // Call auto-tune if 'T' is received
        }
    }
}

// --- Core Functions ---

PhaseAnalysisResult analyzePhases() {
    PhaseAnalysisResult result;
    float a = xk1(0), b = xk1(1), c = xk1(2);

    a -= Config::TQA; b -= Config::TQA; c -= Config::TQA;
    
    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta  = (2.0f / 3.0f) * (Config::SQ3 / 2.0f) * (b - c);
    
    result.angle = atan2f(beta, alpha);
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


// --- Utility and Support Functions (Unchanged) ---

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
            pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW); TCNT1 = 0;
            analogWrite(Config::CP, EPW); delayMicroseconds(Config::EPWU);
            analogWrite(Config::CP, 0); pinMode(Config::CP, INPUT);
            //delayMicroseconds(Config::FSTU);
            rT[0] += cAR(Config::PA); rT[1] += cAR(Config::PB); rT[2] += cAR(Config::PC);
        }
        BLA::Matrix<NKF1, 1> avg_rT;
        for (int i = 0; i < 3; i++) avg_rT(i) = static_cast<float>(rT[i]) / static_cast<float>(Config::DNS);
        aMKF(avg_rT);
        int cAA = static_cast<int>((xk1(0) + xk1(1) + xk1(2)) / 3.0f);
        int err = cAA - Config::TQA;
        if (abs(err) <= Config::PCT) { Serial.println("--- Auto-tune complete! ---"); break; }
        EPW += (err < 0) ? Config::PAS : -Config::PAS;
        EPW = constrain(EPW, 1, 255);
        delay(5);
    }
    Serial.print("Final Tuned PWM Value: "); Serial.println(EPW);
}

int gDNS(float cV) {
    float aV = fabsf(cV);
    if (aV < Config::LV_TH) return Config::MAX_S;
    if (aV > Config::HV_TH) return Config::MIN_S;
    float nV = (aV - Config::LV_TH) / (Config::HV_TH - Config::LV_TH);
    int nS = Config::MAX_S - static_cast<int>(nV * (Config::MAX_S - Config::MIN_S));
    return constrain(nS, Config::MIN_S, Config::MAX_S);
}

void eSC() {
    long r[3] = {0};
    int nST = gDNS(estimatedVelocity);
    for (int i = 0; i < nST; i++) {
        pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW); TCNT1 = 0;
        analogWrite(Config::CP, EPW); delayMicroseconds(Config::EPWU);
        analogWrite(Config::CP, 0); pinMode(Config::CP, INPUT);
        //delayMicroseconds(Config::FSTU);
        r[0] += cAR(Config::PA); r[1] += cAR(Config::PB); r[2] += cAR(Config::PC);
    }
    BLA::Matrix<NKF1, 1> avg_r;
    for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(r[i]) / static_cast<float>(nST);
    aMKF(avg_r);
    pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW);
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

int cAR(int p) {
    return analogRead(p);
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
Thought for 113s
I'll first define a small todo list for adding a non-movement confidence layer and a startup rotation detector, then implement the code changes by creating a new Arduino sketch file with the updated logic and logging.

Updated todo list and thought for 37s
I'll create a short task list and start by implementing the non-movement confidence layer, then add the startup rotation detector, integrate them into the state machine, and update logging.

Updated todo list and used apply_patch and thought for 41s
I'll add a motion confidence layer and startup detector, integrate them into your state machine, and write a new Arduino sketch file containing the changes. Then I’ll mark the tasks done.

Edited rotation_encoder_confidence.ino
Thought for 24s
I’ll add a non-movement confidence layer and a startup rotation detector, then integrate them into the state machine and logging.

Drop-in sketch (with confidence and startup layers)
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
    constexpr int DNS = 8; // DEFAULT_NUM_SAMPLES
    constexpr int MIN_S = 4; // MIN_NUM_SAMPLES
    constexpr int MAX_S = 8; // MAX_NUM_SAMPLES
    constexpr float LV_TH = 5.0f; // LOW_VELOCITY_THRESHOLD (rad/s)
    constexpr float HV_TH = 15.0f; // HIGH_VELOCITY_THRESHOLD (rad/s)

    // Excitation & Auto-Tuning
    constexpr int EPWU = 512, FSTU = 1;
    constexpr int TQA = 32, PAS = 1, PCI = 2000, PCT = 1;

    // Kalman Filter (for ADC smoothing)
    constexpr float QDK1 = 1.0f, RDK1 = 80.1f;
    constexpr float MATRIX_SINGULARITY_THRESHOLD = 1e-5f;

    // --- Core Estimation & Tracking Parameters ---
    constexpr float MIN_SIGNAL_MAGNITUDE = 1.0f; // Minimum vector magnitude for a reliable signal (in ADC units). This is a critical tuning parameter.
    constexpr float VST = 0.01f; // Minimum velocity to consider moving (rad/s).
    constexpr int ENGAGE_CONFIRMATION_FRAMES = 10; // How many frames of consistent motion are needed to lock tracking.

    // --- Motion Confidence Layer ---
    constexpr float OMEGA_NOISE_TH = 0.5f;         // optional noise floor (rad/s) for plane omega (not strictly required)
    constexpr float STARTUP_OMEGA_TH = 2.0f;       // rad/s required to consider rotation startup
    constexpr int STARTUP_CONFIRM_FRAMES = 6;      // frames to confirm startup direction
    constexpr float STATIONARY_ENTER_OMEGA = 0.2f; // rad/s under which we consider stationary accumulation
    constexpr float STATIONARY_EXIT_OMEGA  = 0.4f; // rad/s above which we exit stationary state
    constexpr float CONF_ALPHA = 0.2f;             // smoothing for confidence signals (0..1)
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
float lastRawElectricalAngle = 0.0f;
float unwrappedMeasuredAngle = 0.0f;
float estimatedAngle = 0.0f;
float estimatedVelocity = 0.0f;

// Timestamps and values for velocity calculation
unsigned long pVTs = 0;
unsigned long last_update_timestamp = 0;
float last_unwrapped_angle_for_vel = 0.0f;

// Mechanical Position State
long rC = 0;  // revolution_count
long eMS = 0; // estimated_mechanical_steps

// --- Motion Confidence Layer State ---
float stationaryConfidence = 1.0f;  // 1.0 = very confident stationary, 0.0 = moving
float startupConfidence = 0.0f;     // 0.0..1.0 confidence of rotation startup in a direction
int startupDirection = 0;           // -1, 0, +1 (electrical rotation direction)
int startupConsistentFrames = 0;
float planeOmega = 0.0f;            // instantaneous angular velocity in alpha-beta plane (rad/s)

float lastAlpha = 0.0f, lastBeta = 0.0f;
unsigned long lastAlphaBetaTs = 0;

// --- Structs ---
struct PhaseAnalysisResult {
    float angle = 0.0f;
    float magnitude = 0.0f;
    float alpha = 0.0f;
    float beta = 0.0f;
    bool is_reliable = false;
};

// --- Function Prototypes ---
void initM();
void aTEP();
void eSC();
PhaseAnalysisResult analyzePhases();
void updateUnwrappedMeasuredAngle(float currentRawAngle);
void lD(const PhaseAnalysisResult& result);
int cAR(int);
void aMKF(const BLA::Matrix<NKF1, 1>& z);
int gDNS(float velocity);
unsigned long safeMicrosDiff(unsigned long current, unsigned long previous);
void updateMotionConfidence(const PhaseAnalysisResult& result, unsigned long now);

unsigned long safeMicrosDiff(unsigned long current, unsigned long previous) {
    return (current >= previous) ? (current - previous) : ((ULONG_MAX - previous) + current + 1);
}

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println("Initializing Corrected BLDC Sensorless Encoder (with confidence layers)...");

    initM();
    analogReference(INTERNAL);
    ADCSRA = (ADCSRA & ~0x07) | 0x05;
    TCCR1B = (TCCR1B & 0b11111000) | 0x01;

    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW);
    aTEP();

    // Initialize ADC Kalman Filter
    xk1.Fill(static_cast<float>(Config::TQA));
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 100.0f;
    
    // Initialize timestamps
    pVTs = micros();
    last_update_timestamp = pVTs;
    lastAlphaBetaTs = pVTs;

    Serial.println("System ready. Current state: UNLOCKED");
}

void loop() {
    unsigned long now = micros();

    eSC(); // Perform ADC sampling cycle
    PhaseAnalysisResult reading = analyzePhases(); // Get angle and signal quality
    updateMotionConfidence(reading, now); // Update planeOmega, stationary/startup confidences

    // --- Core State Machine ---
    switch (currentState) {
        case UNLOCKED: {
            estimatedVelocity = 0.0f;
            consistentMovementCounter = 0; // Reset counter in UNLOCKED state

            // Gate engagement on startup detection and non-stationary confidence
            bool can_engage = reading.is_reliable &&
                              (startupConfidence > 0.7f || fabsf(planeOmega) > Config::VST) &&
                              (stationaryConfidence < 0.7f);

            if (can_engage) {
                Serial.println("UNLOCKED -> ENGAGING: Reliable signal and startup detected.");
                
                // Re-synchronize unwrappedMeasuredAngle with the current reading.angle
                float currentWrappedAngle = reading.angle; // The raw angle from -PI to PI
                float lastUnwrappedRemainder = fmodf(unwrappedMeasuredAngle, Config::TPF);
                if (lastUnwrappedRemainder < 0) {
                    lastUnwrappedRemainder += Config::TPF;
                }
                float angleDiffWrapped = currentWrappedAngle - lastUnwrappedRemainder;
                if (angleDiffWrapped > Config::PF) angleDiffWrapped -= Config::TPF;
                else if (angleDiffWrapped < -Config::PF) angleDiffWrapped += Config::TPF;
                unwrappedMeasuredAngle += angleDiffWrapped;
                lastRawElectricalAngle = reading.angle;
                
                last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
                last_update_timestamp = now;

                // Seed velocity with planeOmega sign and magnitude
                estimatedVelocity = planeOmega;

                currentState = ENGAGING;
            }
            break;
        }
        
        case ENGAGING: {
            // Output zero velocity while engaging to avoid jumps on consumers
            estimatedVelocity = 0.0f;

            if (!reading.is_reliable) {
                Serial.println("ENGAGING -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                consistentMovementCounter = 0;
                break;
            }

            // Keep lastRawElectricalAngle aligned with the current raw reading.
            lastRawElectricalAngle = reading.angle;

            // Instantaneous velocity proxy uses planeOmega (unwrap-independent)
            float instant_vel = planeOmega;
            
            // Count consistent motion only when not stationary
            if (fabsf(instant_vel) > Config::VST && stationaryConfidence < 0.7f) {
                consistentMovementCounter++;
            } else {
                consistentMovementCounter = 0;
            }

            last_update_timestamp = now;

            // Promote to LOCKED with either enough consistent frames or strong startup confidence
            if ((consistentMovementCounter > Config::ENGAGE_CONFIRMATION_FRAMES && reading.is_reliable) ||
                startupConfidence > 0.85f) {
                Serial.println("ENGAGING -> LOCKED: Motion confirmed.");
                currentState = LOCKED;
                last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
                estimatedAngle = unwrappedMeasuredAngle;
                estimatedVelocity = instant_vel; // seed
            }
            break;
        }

        case LOCKED: {
            if (!reading.is_reliable) {
                Serial.println("LOCKED -> UNLOCKED: Signal lost.");
                currentState = UNLOCKED;
                estimatedVelocity = 0.0f;
                break;
            }

            // Freeze accumulation near stationary to prevent drift and spurious reversals
            if (stationaryConfidence > 0.8f) {
                estimatedVelocity = 0.0f;
                last_update_timestamp = now;
                lastRawElectricalAngle = reading.angle; // track wrapped input
                estimatedAngle = unwrappedMeasuredAngle; // hold
                break;
            }

            // ONLY update unwrappedMeasuredAngle with the change during the LOCKED state
            updateUnwrappedMeasuredAngle(reading.angle);
            
            // Calculate instantaneous velocity
            float dt_sec = safeMicrosDiff(now, last_update_timestamp) / 1e6f;
            float instantaneous_velocity = 0.0f;
            if (dt_sec > 1e-5f) {
                instantaneous_velocity = (unwrappedMeasuredAngle - last_unwrapped_angle_for_vel) / dt_sec;
            }
            last_update_timestamp = now;
            last_unwrapped_angle_for_vel = unwrappedMeasuredAngle;
            
            estimatedVelocity = instantaneous_velocity;
            estimatedAngle = unwrappedMeasuredAngle;
            break;
        }
    }

    // Update final derived values for output
    eMS = round(estimatedAngle / Config::TPF * Config::ESR);
    rC = round(estimatedAngle / Config::TPF);
    
    pVTs = now;
    lD(reading); // Log data to Serial Plotter

    if (Serial.available() > 0) {
        char cmd = toupper(Serial.read());
        if (cmd == 'T') {
            aTEP(); // Call auto-tune if 'T' is received
        }
    }
}

// --- Core Functions ---

PhaseAnalysisResult analyzePhases() {
    PhaseAnalysisResult result;
    float a = xk1(0), b = xk1(1), c = xk1(2);

    a -= Config::TQA; b -= Config::TQA; c -= Config::TQA;
    
    float alpha = (2.0f / 3.0f) * (a - 0.5f * (b + c));
    float beta  = (2.0f / 3.0f) * (Config::SQ3 / 2.0f) * (b - c);
    
    result.angle = atan2f(beta, alpha);
    result.magnitude = sqrtf(alpha * alpha + beta * beta);
    result.alpha = alpha;
    result.beta = beta;
    result.is_reliable = result.magnitude > Config::MIN_SIGNAL_MAGNITUDE;
    return result;
}

void updateMotionConfidence(const PhaseAnalysisResult& result, unsigned long now) {
    // Compute derivatives in alpha-beta plane and infer angular velocity independent of unwrap
    if (lastAlphaBetaTs == 0) {
        lastAlpha = result.alpha;
        lastBeta = result.beta;
        lastAlphaBetaTs = now;
        planeOmega = 0.0f;
        return;
    }

    float dt_sec = safeMicrosDiff(now, lastAlphaBetaTs) / 1e6f;
    if (dt_sec < 1e-6f) dt_sec = 1e-6f;

    float dalpha = (result.alpha - lastAlpha) / dt_sec;
    float dbeta  = (result.beta  - lastBeta)  / dt_sec;

    // omega = (x*y' - y*x') / (x^2 + y^2)
    float denom = (result.alpha * result.alpha + result.beta * result.beta);
    if (denom < 1e-6f) denom = 1e-6f;
    float rawPlaneOmega = (result.alpha * dbeta - result.beta * dalpha) / denom;

    // One-pole smoothing
    const float omega_alpha = 0.3f;
    planeOmega = omega_alpha * rawPlaneOmega + (1.0f - omega_alpha) * planeOmega;

    // --- Stationary confidence update (with hysteresis) ---
    float stationary_indicator = 0.0f;
    if (fabsf(planeOmega) < Config::STATIONARY_ENTER_OMEGA || result.magnitude < Config::MIN_SIGNAL_MAGNITUDE * 1.25f) {
        stationary_indicator = 1.0f;
    } else if (fabsf(planeOmega) > Config::STATIONARY_EXIT_OMEGA) {
        stationary_indicator = 0.0f;
    } else {
        stationary_indicator = stationaryConfidence;
    }
    stationaryConfidence = Config::CONF_ALPHA * stationary_indicator + (1.0f - Config::CONF_ALPHA) * stationaryConfidence;

    // --- Startup detector update ---
    if (fabsf(planeOmega) > Config::STARTUP_OMEGA_TH && stationaryConfidence < 0.5f) {
        int sign = (planeOmega > 0.0f) ? 1 : -1;
        if (sign == startupDirection || startupDirection == 0) {
            startupConsistentFrames++;
        } else {
            startupConsistentFrames = 1;
        }
        startupDirection = sign;
    } else {
        if (startupConsistentFrames > 0) startupConsistentFrames--;
        if (abs(startupConsistentFrames) == 0) startupDirection = 0;
    }
    float targetStartupConf = (startupConsistentFrames >= Config::STARTUP_CONFIRM_FRAMES)
        ? 1.0f
        : (float)startupConsistentFrames / (float)Config::STARTUP_CONFIRM_FRAMES;
    if (targetStartupConf < 0.0f) targetStartupConf = 0.0f;
    if (targetStartupConf > 1.0f) targetStartupConf = 1.0f;
    startupConfidence = Config::CONF_ALPHA * targetStartupConf + (1.0f - Config::CONF_ALPHA) * startupConfidence;

    lastAlpha = result.alpha;
    lastBeta = result.beta;
    lastAlphaBetaTs = now;
}

void updateUnwrappedMeasuredAngle(float currentRawAngle) {
    float angleDelta = currentRawAngle - lastRawElectricalAngle;
    if (angleDelta > Config::PF) angleDelta -= Config::TPF;
    else if (angleDelta < -Config::PF) angleDelta += Config::TPF;
    
    unwrappedMeasuredAngle += angleDelta;
    lastRawElectricalAngle = currentRawAngle;
}


// --- Utility and Support Functions (Unchanged) ---

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
            pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW); TCNT1 = 0;
            analogWrite(Config::CP, EPW); delayMicroseconds(Config::EPWU);
            analogWrite(Config::CP, 0); pinMode(Config::CP, INPUT);
            //delayMicroseconds(Config::FSTU);
            rT[0] += cAR(Config::PA); rT[1] += cAR(Config::PB); rT[2] += cAR(Config::PC);
        }
        BLA::Matrix<NKF1, 1> avg_rT;
        for (int i = 0; i < 3; i++) avg_rT(i) = static_cast<float>(rT[i]) / static_cast<float>(Config::DNS);
        aMKF(avg_rT);
        int cAA = static_cast<int>((xk1(0) + xk1(1) + xk1(2)) / 3.0f);
        int err = cAA - Config::TQA;
        if (abs(err) <= Config::PCT) { Serial.println("--- Auto-tune complete! ---"); break; }
        EPW += (err < 0) ? Config::PAS : -Config::PAS;
        EPW = constrain(EPW, 1, 255);
        delay(5);
    }
    Serial.print("Final Tuned PWM Value: "); Serial.println(EPW);
}

int gDNS(float cV) {
    float aV = fabsf(cV);
    if (aV < Config::LV_TH) return Config::MAX_S;
    if (aV > Config::HV_TH) return Config::MIN_S;
    float nV = (aV - Config::LV_TH) / (Config::HV_TH - Config::LV_TH);
    int nS = Config::MAX_S - static_cast<int>(nV * (Config::MAX_S - Config::MIN_S));
    return constrain(nS, Config::MIN_S, Config::MAX_S);
}

void eSC() {
    long r[3] = {0};
    int nST = gDNS(estimatedVelocity);
    for (int i = 0; i < nST; i++) {
        pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW); TCNT1 = 0;
        analogWrite(Config::CP, EPW); delayMicroseconds(Config::EPWU);
        analogWrite(Config::CP, 0); pinMode(Config::CP, INPUT);
        //delayMicroseconds(Config::FSTU);
        r[0] += cAR(Config::PA); r[1] += cAR(Config::PB); r[2] += cAR(Config::PC);
    }
    BLA::Matrix<NKF1, 1> avg_r;
    for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(r[i]) / static_cast<float>(nST);
    aMKF(avg_r);
    pinMode(Config::CP, OUTPUT); digitalWrite(Config::CP, LOW);
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

int cAR(int p) {
    return analogRead(p);
}

void lD(const PhaseAnalysisResult& result) {
    // For Arduino IDE Serial Plotter
    Serial.print(xk1(0)-Config::TQA); Serial.print(",");
    Serial.print(xk1(1)-Config::TQA); Serial.print(",");
    Serial.print(xk1(2)-Config::TQA); Serial.print(",");
    Serial.print("Angle:"); Serial.print(estimatedAngle, 2); Serial.print(",");
    Serial.print("Vel:"); Serial.print(estimatedVelocity, 3); Serial.print(",");
    Serial.print("Mag:"); Serial.print(result.magnitude, 2); Serial.print(",");
    Serial.print(eMS); Serial.print(",");
    Serial.print(rC); Serial.print(",");
    Serial.print("pOm:"); Serial.print(planeOmega, 3); Serial.print(",");
    Serial.print("SC:"); Serial.print(stationaryConfidence, 2); Serial.print(",");
    Serial.print("SU:"); Serial.print(startupConfidence, 2); Serial.print(",");
    Serial.print("Dir:"); Serial.print(startupDirection);
    Serial.print("State:"); Serial.println(currentState); // 0=UNLOCKED, 1=ENGAGING, 2=LOCKED
}
