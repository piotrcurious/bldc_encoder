#include <Arduino.h>

// --- Pin Definitions ---
// These pins are connected to the three motor phases (U, V, W).
// They will be dynamically configured as OUTPUT (HIGH/LOW) or INPUT (FLOAT).
const int PHASE_U_PIN = A0; // Analog Pin A0, also digital pin D14
const int PHASE_V_PIN = A1; // Analog Pin A1, also digital pin D15
const int PHASE_W_PIN = A2; // Analog Pin A2, also digital pin D16

// This pin provides the Pulse Width Modulation (PWM) signal for motor speed control.
// Its state (HIGH/LOW) will dictate the polarity of the active phase drive.
const int PWM_SPEED_PIN = 9; // Digital Pin 9 (PWM capable - OC1A on Timer1)

// --- Motor Control Parameters ---
volatile int motorSpeedPWMValue = 100; // PWM duty cycle (0-255) for speed control.
                                       // Adjust this value to change motor RPM.
volatile int commutationStep = 0;      // Current commutation step (0-5, representing the 6 steps)
volatile unsigned long lastCommutationTime = 0; // Timestamp of the last commutation
volatile unsigned long commutationPeriodUs = 2000; // Initial target period in microseconds for commutation.
                                                    // This is adjusted by BEMF sensing.

// Back-EMF (BEMF) threshold for zero-crossing detection.
// For a 5V supply, 512 is approximately VCC/2 (2.5V) on a 10-bit ADC.
// This is the virtual neutral point reference.
const int BEMF_THRESHOLD = 512;

// Number of open-loop commutations during the startup phase.
// This helps the motor gain initial momentum and align before switching to sensorless control.
const int STARTUP_COMMUTATIONS = 30; // Motor will perform 30 full rotations open loop (30 * 6 steps)
volatile int startupCounter = 0;     // Counter for open-loop startup steps.

volatile bool motorStarted = false;  // Flag to indicate if open-loop startup is complete.
volatile bool currentPwmStateHigh = false; // Tracks the state of PWM Pin 9 (HIGH or LOW part of cycle)
volatile bool bemfCheckReady = false; // Flag to indicate when BEMF sensing should occur.

// --- Commutation Table ---
// This structure defines the state of the three motor phases for each of the 6 commutation steps.
// - highPin: The motor phase pin conceptually connected to VCC (or actively driven HIGH).
// - lowPin: The motor phase pin conceptually connected to GND (or actively driven LOW).
// - floatPin: The motor phase pin to set to INPUT (floating) for Back-EMF sensing.
// - sensePinAnalog: The analog input pin corresponding to the floatPin for reading BEMF.
// - bemfExpectedDirection: 1 for rising BEMF, -1 for falling BEMF during zero-crossing.
struct CommutationState {
  int highPin;
  int lowPin;
  int floatPin;
  int sensePinAnalog;
  int bemfExpectedDirection; // 1 for rising, -1 for falling BEMF
};

// The actual commutation sequence table for 6-step trapezoidal commutation.
// This table defines the *forward* drive polarity.
CommutationState commutationTable[] = {
  //  High Phase,     Low Phase,      Float Phase,    Analog Pin, Expected BEMF Direction
  {PHASE_U_PIN,   PHASE_V_PIN,     PHASE_W_PIN,    A2, -1}, // Step 0: U-HIGH, V-LOW, W-FLOAT (Sense W, falling)
  {PHASE_U_PIN,   PHASE_W_PIN,     PHASE_V_PIN,    A1,  1}, // Step 1: U-HIGH, W-LOW, V-FLOAT (Sense V, rising)
  {PHASE_V_PIN,   PHASE_W_PIN,     PHASE_U_PIN,    A0, -1}, // Step 2: V-HIGH, W-LOW, U-FLOAT (Sense U, falling)
  {PHASE_V_PIN,   PHASE_U_PIN,     PHASE_W_PIN,    A2,  1}, // Step 3: V-HIGH, U-LOW, W-FLOAT (Sense W, rising)
  {PHASE_W_PIN,   PHASE_U_PIN,     PHASE_V_PIN,    A1, -1}, // Step 4: W-HIGH, U-LOW, V-FLOAT (Sense V, falling)
  {PHASE_W_PIN,   PHASE_V_PIN,     PHASE_U_PIN,    A0,  1}  // Step 5: W-HIGH, V-LOW, U-FLOAT (Sense U, rising)
};

// --- Function Prototypes ---
void setMotorPhase(int highPin, int lowPin, int floatPin, bool applyReversePolarity);
void advanceCommutationStep();
void handleBEMFSensing(); // Now called from ISR

// --- Interrupt Service Routines (ISRs) ---

// This ISR fires when Timer1 counter (TCNT1) overflows (reaches TOP, then resets to 0).
// In Fast PWM mode (Mode 14), this happens at the start of a new PWM cycle,
// and OC1A (Pin 9) is set HIGH at this point (assuming non-inverting mode).
ISR(TIMER1_OVF_vect) {
  // This marks the beginning of the PWM "ON" time (Pin 9 is HIGH).
  currentPwmStateHigh = true;
  CommutationState currentStep = commutationTable[commutationStep];
  // Apply the *forward* drive for the current commutation step.
  setMotorPhase(currentStep.highPin, currentStep.lowPin, currentStep.floatPin, false);
}

// This ISR fires when Timer1 counter (TCNT1) matches OCR1A.
// In Fast PWM mode (Mode 14), this causes OC1A (Pin 9) to go LOW.
// This marks the beginning of the PWM "OFF" time.
ISR(TIMER1_COMPA_vect) {
  // This marks the beginning of the PWM "OFF" time (Pin 9 is LOW).
  currentPwmStateHigh = false;
  CommutationState currentStep = commutationTable[commutationStep];

  // Apply the *reverse* drive as per the "novel" bidirectional requirement.
  // The previously HIGH pin goes LOW, and the previously LOW pin goes HIGH.
  setMotorPhase(currentStep.highPin, currentStep.lowPin, currentStep.floatPin, true);

  // Set flag to perform BEMF sensing. It's safe to do short, non-blocking
  // analog reads and logic here, as we are in an ISR and need tight timing.
  handleBEMFSensing();
}


void setup() {
  Serial.begin(9600); // Initialize serial communication for debugging
  Serial.println("BLDC Controller Starting with PWM Sync...");

  // Initialize all phase pins as inputs to ensure a safe default state
  pinMode(PHASE_U_PIN, INPUT);
  pinMode(PHASE_V_PIN, INPUT);
  pinMode(PHASE_W_PIN, INPUT);

  // Configure Timer1 for Fast PWM mode (Mode 14 - WGM13, WGM12, WGM11)
  // Pin 9 (OC1A) is used for PWM.
  // Clear OC1A on Compare Match, set on TOP (non-inverting mode).
  TCCR1A = _BV(COM1A1) | _BV(WGM11);
  TCCR1B = _BV(WGM13) | _BV(WGM12); // WGM13 and WGM12 for Mode 14
  
  // Set Timer1 prescaler to 8 (CS11).
  // With 16MHz clock, this results in a PWM frequency of (16MHz / 8) / (ICR1 + 1).
  // If ICR1 is 255 (8-bit like analogWrite), freq = 2MHz / 256 = ~7.8 kHz.
  TCCR1B |= _BV(CS11);

  // Set TOP value for PWM. Using 255 for 8-bit resolution, compatible with analogWrite.
  ICR1 = 255; // TOP value for Fast PWM.

  // Set the initial PWM duty cycle. This value will be loaded into OCR1A.
  OCR1A = motorSpeedPWMValue;

  // Enable Timer1 Compare Match A Interrupt (fires when TCNT1 == OCR1A)
  TIMSK1 |= _BV(OCIE1A);
  // Enable Timer1 Overflow Interrupt (fires when TCNT1 reaches TOP, then resets to 0)
  TIMSK1 |= _BV(TOIE1);

  // Set the PWM pin as an output.
  pinMode(PWM_SPEED_PIN, OUTPUT);

  // Initial setup of phase pins (all float)
  setMotorPhase(0, 0, 0, false); // All floating initially.

  // The main loop will now largely be idle, as commutation and BEMF sensing
  // are driven by Timer1 interrupts.
}

void loop() {
  // The main loop is now primarily for debugging output or handling higher-level tasks
  // not requiring precise timing (e.g., serial command input to change speed).

  // Example: Change speed via serial input
  if (Serial.available()) {
    int newSpeed = Serial.parseInt();
    if (newSpeed >= 0 && newSpeed <= 255) {
      noInterrupts(); // Temporarily disable interrupts for shared variable access
      motorSpeedPWMValue = newSpeed;
      OCR1A = motorSpeedPWMValue; // Update PWM duty cycle
      interrupts(); // Re-enable interrupts
      Serial.print("Motor Speed set to: ");
      Serial.println(motorSpeedPWMValue);
    }
  }

  // Small delay to prevent continuous printing, if needed for other debug.
  // delay(10);
}

/**
 * @brief Configures the states (HIGH, LOW, FLOAT) for the three motor phase pins.
 *
 * This function is called from within the Timer1 ISRs to dynamically set up the pins
 * based on the current commutation step and the state of the PWM signal (Pin 9).
 *
 * @param highPin The digital pin number for the phase to be set HIGH.
 * @param lowPin The digital pin number for the phase to be set LOW.
 * @param floatPin The digital pin number for the phase to be set to INPUT (floating).
 * @param applyReversePolarity If true, the highPin will be set LOW and lowPin will be set HIGH.
 * This implements the "novel" bidirectional current flow during PWM OFF.
 */
void setMotorPhase(int highPin, int lowPin, int floatPin, bool applyReversePolarity) {
  // Reset all phase pins to INPUT mode first to safely transition states
  // and prevent momentary short circuits or unintended current paths.
  pinMode(PHASE_U_PIN, INPUT);
  pinMode(PHASE_V_PIN, INPUT);
  pinMode(PHASE_W_PIN, INPUT);

  if (highPin != 0 && lowPin != 0 && floatPin != 0) { // Only apply if valid pins are provided
    if (!applyReversePolarity) {
      // Normal forward drive: highPin to VCC, lowPin to GND.
      pinMode(highPin, OUTPUT);
      digitalWrite(highPin, HIGH);

      pinMode(lowPin, OUTPUT);
      digitalWrite(lowPin, LOW);
    } else {
      // "Novel" reverse drive (during PWM OFF): highPin to GND, lowPin to VCC.
      // This aims to actively reverse current or provide active freewheeling.
      pinMode(highPin, OUTPUT);
      digitalWrite(highPin, LOW);

      pinMode(lowPin, OUTPUT);
      digitalWrite(lowPin, HIGH);
    }
    // The float phase remains floating for BEMF sensing.
    pinMode(floatPin, INPUT);
  }
}


/**
 * @brief Advances the motor to the next commutation step in the sequence.
 *
 * This function updates the `commutationStep` counter (circularly from 0 to 5).
 * It also prints the current step to the serial monitor for debugging.
 * This is now called from within the `handleBEMFSensing` function, which is ISR-driven.
 */
void advanceCommutationStep() {
  commutationStep = (commutationStep + 1) % 6; // Move to the next step (0-5 cycle)
  lastCommutationTime = micros(); // Record the time of this commutation

  // For debugging, consider sending this via Serial only periodically or via a flag
  // as Serial.print is slow and can block ISRs.
  // Serial.print("Commutating to step: ");
  // Serial.println(commutationStep);
}

/**
 * @brief Performs Back-EMF (BEMF) sensing on the floating motor phase.
 *
 * This function is called from the TIMER1_COMPA_vect ISR (when PWM goes LOW).
 * It reads the analog value from the floating phase and checks for a zero-crossing.
 *
 * If a zero-crossing is detected (indicating the optimal time to commutate),
 * the motor is advanced to the next step, and the commutation period is adjusted
 * for closed-loop speed control.
 *
 * NOTE: BEMF sensing accuracy might be compromised due to the active "reverse drive"
 * during the PWM OFF time. This implementation assumes the BEMF is still discernible.
 */
void handleBEMFSensing() {
  CommutationState currentStep = commutationTable[commutationStep];
  int bemfValue = analogRead(currentStep.sensePinAnalog);
  bool zeroCrossingDetected = false;

  // --- Startup Sequence (Open Loop Control) ---
  if (!motorStarted) {
    // In startup, we advance based on a fixed period (commutationPeriodUs).
    if (micros() - lastCommutationTime >= commutationPeriodUs) {
      advanceCommutationStep();
      startupCounter++;
      // Gradually decrease the period to accelerate the motor.
      if (commutationPeriodUs > 500) { // Keep a minimum period
        commutationPeriodUs -= 10;
      }
      if (startupCounter >= STARTUP_COMMUTATIONS * 6) {
        motorStarted = true; // Startup complete, switch to closed-loop
        Serial.println("Startup complete. Switching to sensorless control.");
      }
    }
    return; // Don't do BEMF sensing during open-loop startup.
  }

  // --- Closed Loop Operation (Sensorless BEMF Control) ---
  // The BEMF detection is performed when Pin 9 is LOW (active reverse drive).
  // The expected direction of zero-crossing must be checked against the
  // virtual neutral point (BEMF_THRESHOLD).

  // Serial.print("Sensing BEMF on A"); Serial.print(currentStep.sensePinAnalog - A0);
  // Serial.print(": "); Serial.println(bemfValue);

  // Simplified zero-crossing detection.
  // This logic must consider the BEMF_THRESHOLD relative to the expected
  // BEMF waveform polarity for the *floating* phase during the active reverse drive.
  // If BEMF_THRESHOLD is mid-supply, and BEMF should fall/rise towards it:
  if (currentStep.bemfExpectedDirection == -1) { // Expecting falling BEMF
    if (bemfValue < BEMF_THRESHOLD) {
      zeroCrossingDetected = true;
    }
  } else if (currentStep.bemfExpectedDirection == 1) { // Expecting rising BEMF
    if (bemfValue > BEMF_THRESHOLD) {
      zeroCrossingDetected = true;
    }
  }

  if (zeroCrossingDetected) {
    unsigned long currentTime = micros();
    unsigned long timeElapsed = currentTime - lastCommutationTime;

    // Adjust the commutation period based on sensed time.
    // Simple proportional control: if timeElapsed is shorter, motor is faster,
    // so next period should be shorter. If longer, motor is slower, period should be longer.
    // This is a simple smoothing filter for the commutation period.
    commutationPeriodUs = (commutationPeriodUs * 7 + timeElapsed) / 8;

    advanceCommutationStep(); // Commutate to the next step
  }
}
