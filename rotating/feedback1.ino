#include <Arduino.h>

// --- Pin Definitions ---
// These pins are connected to the three motor phases (U, V, W).
// They will be dynamically configured as OUTPUT (HIGH/LOW) or INPUT (FLOAT).
const int PHASE_U_PIN = A0; // Analog Pin A0, also digital pin D14
const int PHASE_V_PIN = A1; // Analog Pin A1, also digital pin D15
const int PHASE_W_PIN = A2; // Analog Pin A2, also digital pin D16

// This pin provides the Pulse Width Modulation (PWM) signal for motor speed control.
// It's conceptually tied into the low-side drive of the active phase in each commutation step.
const int PWM_SPEED_PIN = 9; // Digital Pin 9 (PWM capable)

// --- Motor Control Parameters ---
int motorSpeed = 100; // Initial PWM value (0-255) for speed control.
                      // Adjust this value to change motor RPM.
volatile int commutationStep = 0; // Current commutation step (0-5, representing the 6 steps)
unsigned long lastCommutationTime = 0; // Timestamp of the last commutation
unsigned long commutationDelay = 2000; // Initial delay in microseconds for open-loop startup.
                                       // This value will decrease as the motor speeds up.

// Back-EMF (BEMF) threshold for zero-crossing detection.
// For a 5V supply, 512 is approximately VCC/2 (2.5V) on a 10-bit ADC.
// This is the virtual neutral point reference.
const int BEMF_THRESHOLD = 512;

// Number of open-loop commutations during the startup phase.
// This helps the motor gain initial momentum and align before switching to sensorless control.
const int STARTUP_COMMUTATIONS = 30; // Motor will perform 30 full rotations open loop (30 * 6 steps)

// --- Commutation Table ---
// This structure defines the state of the three motor phases for each of the 6 commutation steps.
// - highPin: The motor phase pin to set HIGH (connected to VCC).
// - lowPin: The motor phase pin to set LOW. This pin's connection to ground is conceptually
//           modulated by the PWM_SPEED_PIN (Pin 9) for speed control.
// - floatPin: The motor phase pin to set to INPUT (floating) for Back-EMF sensing.
// - sensePinAnalog: The analog input pin corresponding to the floatPin for reading BEMF.
struct CommutationState {
  int highPin;
  int lowPin;
  int floatPin;
  int sensePinAnalog;
};

// The actual commutation sequence table.
// This implements a standard 6-step trapezoidal commutation pattern.
CommutationState commutationTable[] = {
  //  High Phase,     Low Phase (PWM'd),  Float Phase,    Analog Pin to Sense
  {PHASE_U_PIN,   PHASE_V_PIN,     PHASE_W_PIN,    A2}, // Step 0: U-HIGH, V-LOW (PWM), W-FLOAT (Sense W)
  {PHASE_U_PIN,   PHASE_W_PIN,     PHASE_V_PIN,    A1}, // Step 1: U-HIGH, W-LOW (PWM), V-FLOAT (Sense V)
  {PHASE_V_PIN,   PHASE_W_PIN,     PHASE_U_PIN,    A0}, // Step 2: V-HIGH, W-LOW (PWM), U-FLOAT (Sense U)
  {PHASE_V_PIN,   PHASE_U_PIN,     PHASE_W_PIN,    A2}, // Step 3: V-HIGH, U-LOW (PWM), W-FLOAT (Sense W)
  {PHASE_W_PIN,   PHASE_U_PIN,     PHASE_V_PIN,    A1}, // Step 4: W-HIGH, U-LOW (PWM), V-FLOAT (Sense V)
  {PHASE_W_PIN,   PHASE_V_PIN,     PHASE_U_PIN,    A0}  // Step 5: W-HIGH, V-LOW (PWM), U-FLOAT (Sense U)
};

// --- Function Prototypes ---
void setMotorPhase(int highPin, int lowPin, int floatPin);
void advanceCommutationStep();
void performBEMFSensing();
void setMotorSpeed(int speed);

void setup() {
  Serial.begin(9600); // Initialize serial communication for debugging
  Serial.println("BLDC Controller Starting...");

  // Initialize all phase pins as inputs to ensure a safe default state
  pinMode(PHASE_U_PIN, INPUT);
  pinMode(PHASE_V_PIN, INPUT);
  pinMode(PHASE_W_PIN, INPUT);

  // Set the PWM pin (Pin 9) as an output. This pin will continuously output the PWM signal.
  pinMode(PWM_SPEED_PIN, OUTPUT);

  // Apply the initial motor speed PWM. This will affect the "LOW" phases.
  analogWrite(PWM_SPEED_PIN, motorSpeed);

  // Note on PWM Synchronization:
  // For highly precise BEMF sensing synchronized with the PWM signal,
  // one would typically configure a Timer interrupt (e.g., Timer1, which drives Pin 9 and 10)
  // to trigger the `performBEMFSensing()` function at a specific point in the PWM cycle
  // (ideally, when the PWM signal is LOW, i.e., during the "off" time).
  // Example for Timer1 Overflow Interrupt (fires when TCNT1 wraps around from TOP to 0):
  // TCCR1A = _BV(COM1A1) | _BV(WGM11); // Fast PWM mode, non-inverting output on OC1A (Pin 9)
  // TCCR1B = _BV(WGM13) | _BV(WGM12) | _BV(CS11); // Prescaler 8 (for 31kHz PWM @ 16MHz CPU)
  // ICR1 = 255; // Set TOP to 255 for 8-bit PWM (compatible with analogWrite)
  // TIMSK1 |= _BV(TOIE1); // Enable Timer1 Overflow Interrupt
  // However, for this "mini" and "novel" simplified concept, we will use a small delay
  // in `performBEMFSensing()` to approximate this, making the code simpler.
}

void loop() {
  // --- Startup Sequence (Open Loop Control) ---
  // During startup, the motor is commutated based on a fixed delay,
  // gradually decreasing the delay to accelerate the motor.
  if (commutationStep < STARTUP_COMMUTATIONS * 6) { // Check if still in startup phase
    if (micros() - lastCommutationTime >= commutationDelay) {
      advanceCommutationStep(); // Move to the next commutation step
      lastCommutationTime = micros(); // Record the time of this commutation

      // Gradually decrease the commutation delay to accelerate the motor
      if (commutationDelay > 500) { // Keep a minimum delay to prevent excessively high speed
        commutationDelay -= 10; // Reduce delay by a small amount each step
      }
    }
  }
  // --- Closed Loop Operation (Sensorless BEMF Control) ---
  // Once the startup sequence is complete, the motor operates in closed-loop mode.
  // Commutations are now triggered by detecting Back-EMF zero-crossings.
  else {
    performBEMFSensing(); // Continuously check for BEMF zero-crossing
  }
}

/**
 * @brief Configures the states (HIGH, LOW, FLOAT) for the three motor phase pins.
 *
 * This function is called by `advanceCommutationStep()` to set up the pins
 * for the current commutation step. It ensures that the previous pin states
 * are safely reset before applying the new configuration.
 *
 * @param highPin The digital pin number for the phase to be set HIGH.
 * @param lowPin The digital pin number for the phase to be set LOW. This pin's ground
 * connection is conceptually modulated by the PWM_SPEED_PIN (Pin 9).
 * @param floatPin The digital pin number for the phase to be set to INPUT (floating).
 */
void setMotorPhase(int highPin, int lowPin, int floatPin) {
  // Reset all phase pins to INPUT mode first to safely transition states
  // and prevent momentary short circuits or unintended current paths.
  pinMode(PHASE_U_PIN, INPUT);
  pinMode(PHASE_V_PIN, INPUT);
  pinMode(PHASE_W_PIN, INPUT);

  // 1. Set the HIGH phase: Connects the motor winding to VCC.
  pinMode(highPin, OUTPUT);
  digitalWrite(highPin, HIGH);

  // 2. Set the LOW phase: Connects the motor winding to GND, but conceptually
  //    this path's effectiveness is controlled by the PWM from PWM_SPEED_PIN (Pin 9).
  //    For this "novel" direct Arduino pin approach, we just set the pin LOW.
  pinMode(lowPin, OUTPUT);
  digitalWrite(lowPin, LOW);

  // 3. Set the FLOAT phase: Disconnects the motor winding to allow BEMF to be sensed.
  //    The pin is configured as an INPUT to allow analog reading.
  pinMode(floatPin, INPUT);
}

/**
 * @brief Advances the motor to the next commutation step in the sequence.
 *
 * This function updates the `commutationStep` counter (circularly from 0 to 5)
 * and then calls `setMotorPhase()` to apply the new pin configurations.
 * It also prints the current step to the serial monitor for debugging.
 */
void advanceCommutationStep() {
  commutationStep = (commutationStep + 1) % 6; // Move to the next step (0-5 cycle)
  CommutationState currentStep = commutationTable[commutationStep]; // Get states for the new step

  // Apply the new phase configurations
  setMotorPhase(currentStep.highPin, currentStep.lowPin, currentStep.floatPin);

  // Re-apply the PWM value to the PWM_SPEED_PIN. This ensures the PWM is
  // continuously active and modulating the "LOW" phase.
  analogWrite(PWM_SPEED_PIN, motorSpeed);

  Serial.print("Commutating to step: ");
  Serial.println(commutationStep);
}

/**
 * @brief Performs Back-EMF (BEMF) sensing on the floating motor phase.
 *
 * This is the core of the sensorless control. The function waits for a brief
 * period to allow the BEMF signal to stabilize (and conceptually, for the
 * PWM signal to be in its "off" cycle). It then reads the analog value from
 * the floating phase and checks for a zero-crossing against the BEMF_THRESHOLD.
 *
 * If a zero-crossing is detected (indicating the optimal time to commutate),
 * the motor is advanced to the next step, and the commutation delay is adjusted
 * for closed-loop speed control.
 */
void performBEMFSensing() {
  CommutationState currentStep = commutationTable[commutationStep]; // Get current step info
  int bemfValue = 0;
  bool zeroCrossingDetected = false;

  // Small delay to allow the BEMF voltage to settle after the commutation change,
  // and to ensure we are reading when the PWM signal on Pin 9 is conceptually "LOW"
  // (during its off-time), reducing noise from driving current.
  // This value is critical and depends on motor inductance and PWM frequency.
  delayMicroseconds(20); // Adjust this delay for your specific motor and setup.

  // Read the analog value from the floating phase
  bemfValue = analogRead(currentStep.sensePinAnalog);

  // Uncomment for detailed BEMF debugging:
  // Serial.print("Step "); Serial.print(commutationStep);
  // Serial.print(", BEMF on A"); Serial.print(currentStep.sensePinAnalog - A0);
  // Serial.print(": "); Serial.println(bemfValue);

  // --- Zero-Crossing Detection Logic ---
  // This logic checks for the BEMF crossing the virtual neutral point (BEMF_THRESHOLD)
  // in the expected direction for the current commutation step.
  // This is a simplified zero-crossing detection. A robust solution might involve:
  // - Filtering the BEMF signal (e.g., moving average).
  // - Detecting rising/falling edges across the threshold more explicitly
  //   by tracking the previous BEMF value.
  switch (commutationStep) {
    case 0: // Sense W (A2) after U-HIGH, V-LOW. Expect W falling.
      if (bemfValue < BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
    case 1: // Sense V (A1) after U-HIGH, W-LOW. Expect V rising.
      if (bemfValue > BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
    case 2: // Sense U (A0) after V-HIGH, W-LOW. Expect U falling.
      if (bemfValue < BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
    case 3: // Sense W (A2) after V-HIGH, U-LOW. Expect W rising.
      if (bemfValue > BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
    case 4: // Sense V (A1) after W-HIGH, U-LOW. Expect V falling.
      if (bemfValue < BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
    case 5: // Sense U (A0) after W-HIGH, V-LOW. Expect U rising.
      if (bemfValue > BEMF_THRESHOLD) zeroCrossingDetected = true;
      break;
  }

  if (zeroCrossingDetected) {
    unsigned long currentTime = micros();
    unsigned long timeElapsed = currentTime - lastCommutationTime;

    // Adjust the commutation delay for closed-loop speed control.
    // A simple proportional control: if timeElapsed is shorter, motor is faster,
    // so next delay should be shorter. If longer, motor is slower, delay should be longer.
    // A simple moving average filter helps stabilize the speed.
    commutationDelay = (commutationDelay * 7 + timeElapsed) / 8; // Adjust smoothing factor (7/8) as needed

    advanceCommutationStep(); // Commutate to the next step
    lastCommutationTime = currentTime; // Reset the last commutation time
  }
}

/**
 * @brief Sets the motor speed by adjusting the PWM duty cycle.
 *
 * This function can be called dynamically (e.g., from serial input, a potentiometer)
 * to change the motor's operating speed.
 *
 * @param speed The desired PWM duty cycle (0-255). Values outside this range will be constrained.
 */
void setMotorSpeed(int speed) {
  motorSpeed = constrain(speed, 0, 255); // Constrain speed to valid PWM range
  analogWrite(PWM_SPEED_PIN, motorSpeed); // Apply the new PWM duty cycle
  Serial.print("Motor Speed set to: ");
  Serial.println(motorSpeed);
}

// You could add more functions here, for example:
// - A function to handle serial input to change motor speed.
// - Error handling for stalled motor detection.
// - PID control for more precise speed regulation.
