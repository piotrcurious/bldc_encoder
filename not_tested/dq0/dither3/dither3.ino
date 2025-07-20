 // --- Pin and System Constants ---
const int PHASE_A_PIN = A0;
const int PHASE_B_PIN = A1;
const int PHASE_C_PIN = A2;
const int DITHER_PIN = 9;

const float VREF = 5.0;
const int ADC_MAX = 1023;
const float OFFSET_VOLTAGE = VREF / 2.0;
const float DEG_PER_RAD = 180.0 / PI;

// --- Tuning Parameters ---
const float POS_FILTER_LAMBDA = 0.15;      // Position filter smoothing factor (lower is smoother)
const float VEL_FILTER_ALPHA = 0.1;        // Velocity filter smoothing factor (lower is smoother)
const float MAX_DITHER_AMP = 12.0;         // Max dither amplitude in PWM units
const float DITHER_FREQ = 100.0;           // Dither signal frequency in Hz
const float STARTUP_MAG_THRESHOLD = 0.03;  // Filtered magnitude threshold to start motion
const float STOP_MAG_THRESHOLD = 0.02;     // Filtered magnitude threshold to stop motion (hysteresis)

// --- Data Structures ---
struct Complex {
  float re, im;

  Complex operator+(const Complex& b) const { return {re + b.re, im + b.im}; }
  Complex operator*(float s) const { return {re * s, im * s}; }
  float magnitude() const { return sqrt(re * re + im * im); }
  float angle() const { return atan2(im, re); }
};

// --- Global State Variables ---
Complex z_curr = {0, 0};      // Current raw measurement vector
Complex momentum = {0, 0};    // Filtered position vector

float prev_theta = 0;         // Previous angle for velocity calculation
float velocity = 0;           // Raw velocity (rad/s)
float filtered_velocity = 0;  // Smoothed velocity (rad/s)

unsigned long prev_time = 0;  // For calculating dynamic dt
bool is_running = false;      // State flag for motion

void setup() {
  Serial.begin(115200);
  pinMode(DITHER_PIN, OUTPUT);
  analogWrite(DITHER_PIN, 127); // Set dither PWM to midpoint
  prev_time = millis();
}

void loop() {
  // --- Non-Blocking Timer and Dynamic dt ---
  // Replaces delay() for a more accurate and responsive loop.
  unsigned long current_time = millis();
  float dt = (current_time - prev_time) / 1000.0;
  if (dt <= 0.001) return; // Ensure a minimum time step to avoid division by zero
  prev_time = current_time;

  // --- Adaptive Oversampling ---
  // Decision is now based on the STABLE filtered magnitude, not a noisy raw one.
  float m_mag = momentum.magnitude();
  int oversample_count = 1;
  if (m_mag < 0.05)      oversample_count = 8;
  else if (m_mag < 0.1)  oversample_count = 6;
  else if (m_mag < 0.2)  oversample_count = 4;
  else if (m_mag < 0.5)  oversample_count = 2;
  
  float Va = 0, Vb = 0, Vc = 0;
  for (int i = 0; i < oversample_count; ++i) {
    Va += analogRead(PHASE_A_PIN);
    Vb += analogRead(PHASE_B_PIN);
    Vc += analogRead(PHASE_C_PIN);
  }
  Va = (Va / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  Vb = (Vb / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;
  Vc = (Vc / oversample_count) * VREF / ADC_MAX - OFFSET_VOLTAGE;

  // --- Clarke Transform ---
  float Valpha = (2.0 / 3.0) * (Va - 0.5 * Vb - 0.5 * Vc);
  float Vbeta  = (2.0 / 3.0) * ((sqrt(3.0) / 2.0) * (Vb - Vc));
  z_curr = {Valpha, Vbeta};

  // --- Complex Filtered Momentum (Position Filter) ---
  momentum = momentum * (1.0 - POS_FILTER_LAMBDA) + z_curr * POS_FILTER_LAMBDA;

  // --- State Machine for Startup/Stop ---
  // This logic prevents velocity spikes on startup and zeroes it when stopped.
  if (!is_running && m_mag > STARTUP_MAG_THRESHOLD) {
    is_running = true;
    prev_theta = momentum.angle(); // Initialize angle HERE to prevent jump
  } else if (is_running && m_mag < STOP_MAG_THRESHOLD) {
    is_running = false;
  }

  float theta;
  if (is_running) {
    // --- Velocity Calculation (only when running) ---
    theta = momentum.angle();
    float dtheta = theta - prev_theta;

    // Handle angle wrap-around
    if (dtheta > PI) dtheta -= 2 * PI;
    if (dtheta < -PI) dtheta += 2 * PI;

    velocity = dtheta / dt; // Use dynamic dt
    prev_theta = theta;
  } else {
    // --- System is Stopped ---
    theta = prev_theta; // Hold angle
    velocity = 0;
  }

  // --- Velocity Low-Pass Filter ---
  // A second filter specifically for velocity smooths out noise from differentiation.
  filtered_velocity = filtered_velocity * (1.0 - VEL_FILTER_ALPHA) + velocity * VEL_FILTER_ALPHA;
  if (!is_running) {
    filtered_velocity = 0; // Clamp to zero when stopped
  }
  
  // --- Adaptive Dither ---
  // Amplitude is based on the stable filtered magnitude for robustness.
  float dither_amp = MAX_DITHER_AMP / (1.0 + 50.0 * m_mag);
  float triangle_wave = 2.0 * fabs(fmod(current_time / 1000.0 * DITHER_FREQ, 1.0) - 0.5);
  int pwm_value = 127 + int((triangle_wave - 0.5) * 2.0 * dither_amp);
  analogWrite(DITHER_PIN, pwm_value);

  // --- Output ---
  Serial.print("θ(deg):"); Serial.print(theta * DEG_PER_RAD, 2);
  Serial.print("\tω(rad/s):"); Serial.print(filtered_velocity, 2);
  Serial.print("\t|m|:"); Serial.print(m_mag, 3);
  Serial.print("\tState:"); Serial.print(is_running ? "Run" : "Stop");
  Serial.print("\tOS:"); Serial.println(oversample_count);
}
