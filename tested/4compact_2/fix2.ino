#include <Arduino.h>
#include <BasicLinearAlgebra.h>

// Configuration constants for the sensor system
namespace Config {
    constexpr int PA = A0, PB = A1, PC = A2, CP = 9;
    constexpr float ESR = 1.0f; // ELECTRICAL_STEPS_PER_REV
    const float PF = 3.14159265358979323846f; // PI_F
    constexpr float TPF = 2.0f * PF; // TWO_PI_F
    constexpr float SQ3 = 1.7320508f; // SQRT_3
    
    constexpr int DNS = 16; // DEFAULT_NUM_SAMPLES
    constexpr int MIN_S = 16; // MIN_NUM_SAMPLES
    constexpr int MAX_S = 16; // MAX_NUM_SAMPLES
    constexpr float LV_TH = 5.0f; // LOW_VELOCITY_THRESHOLD
    constexpr float HV_TH = 15.0f; // HIGH_VELOCITY_THRESHOLD

    constexpr float AUD = 0.0002f, VST = 0.0001f; // ANGLE_UPDATE_DEADBAND, VELOCITY_STEP_THRESHOLD
    constexpr unsigned long ZVLT = 200000UL; // ZERO_VEL_LOCK_TIMEOUT
    // BUG FIX: Moved EPW from Config namespace to global scope as it's a state variable, not a config constant.
    constexpr int EPWU = 32, FSTU = 1; // EXCITATION_PULSE_WIDTH_US, FLOATING_SETTLING_TIME_US
    constexpr int TQA = 128, PAS = 1, PCI = 2000, PCT = 1; // TARGET_QUIESCENT_AVG_ADC, PWM_ADJUST_STEP, PWM_CALIBRATION_ITERATIONS, PWM_CALIBRATION_TOLERENCE
    constexpr float QDK1 = 5.0f, RDK1 = 80.0f; // Q_DIAGONAL_VALUE_KF1, R_DIAGONAL_VALUE_KF1
    constexpr float QAE = 1.1f, QVE = 80.1f, RAE = 1.1f; // Q_ANGLE_EKF, Q_VEL_EKF, R_ANGLE_EKF
    
    constexpr int VHS = 4; // VELOCITY_HISTORY_SIZE
    constexpr int AHS = 4; // ANGLE_HISTORY_SIZE
    constexpr int PHS = 4; // PHASE_HISTORY_SIZE
    constexpr float MVC = 0.0003f; // MIN_VELOCITY_CONSISTENCY
    constexpr float VNT = 0.006f; // VELOCITY_NOISE_THRESHOLD
    constexpr float AHNT = 0.0005f; // ANGLE_HISTORY_NOISE_THRESHOLD
    constexpr unsigned long MST = 500UL; // MIN_STABLE_TIME_US
}

int EPW = 1; // EXCITATION_PWM_VALUE (Global state variable)

#define NKF1 3
BLA::Matrix<NKF1, 1> xk1; // x_hat_kf1
BLA::Matrix<NKF1, NKF1> Pk1, Qk1, Rk1, Ik1; // P_kf1, Q_kf1, R_kf1, I_matrix_kf1

#define NEKF 2
BLA::Matrix<NEKF, 1> xe; // x_ekf
BLA::Matrix<NEKF, NEKF> Pe, Qe, Ie; // P_ekf, Q_ekf, I_matrix_ekf
BLA::Matrix<1, 1> Re; // R_ekf

struct HEnt { // HistoryEntry
    float angle;
    float velocity;
    unsigned long timestamp;
};

struct PHEnt { // PhaseHistoryEntry
    float phaseA, phaseB, phaseC;
    unsigned long timestamp;
};

HEnt vH[Config::VHS]; // velocityHistory
HEnt aH[Config::AHS]; // angleHistory
PHEnt pH[Config::PHS]; // phaseHistory
int vHI = 0, aHI = 0, pHI = 0; // velocityHistoryIndex, angleHistoryIndex, phaseHistoryIndex
bool vHF = false, aHF = false, pHF = false; // velocityHistoryFull, angleHistoryFull, phaseHistoryFull

float pREA = 0.0f; // previous_raw_electrical_angle
bool fAR = true; // first_angle_reading

float uEA = 0.0f; // unwrapped_electrical_angle
long rC = 0; // revolution_count
long eMS = 0; // estimated_mechanical_steps

bool eL = false; // ekfLocked
unsigned long lNZM = 0, pEM = 0; // lastNonZeroVelMicros, prev_ekf_micros

float sV = 0.0f; // smoothed_velocity
float vT = 0.0f; // velocity_trend
int cDC = 0; // consistent_direction_count
bool iS = false; // is_stopping
bool iCP = false; // is_consistent_phases

// Function Prototypes
void initM();
void aTEP();
void eSC();
float gIEA();
float nA(float);
void aEKF(float);
void lD();
int cAR(int);
void aMKF(const BLA::Matrix<NKF1, 1>&); // BUG FIX: Changed signature to accept float matrix
void uVH(float, float, unsigned long);
void uAH(float, float, unsigned long);
void uPH(float, float, float, unsigned long);
bool cPC();
bool iVC();
bool iRS();
float gSV();
float gVT();

void uUEA(float);
int gDNS(float);

void kS(float nV, float &sum, float &comp) { // kahanSum
    float y = nV - comp;
    float t = sum + y;
    comp = (t - sum) - y;
    sum = t;
}

void setup() {
    Serial.begin(115200);
    while (!Serial) delay(10);
    Serial.println("Initializing Enhanced BLDC Sensorless Encoder...");
    initM();
    analogReference(INTERNAL);
    ADCSRA = (ADCSRA & ~0x07) | 0x05;
    TCCR1B = (TCCR1B & 0b11111000) | 0x01;
   
    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW);
    aTEP();

    xk1.Fill((float)Config::TQA);
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 100.0f;

    xe.Fill(0.0f);
    Pe.Fill(0.0f);
    for (int i = 0; i < NEKF; i++) Pe(i, i) = 1.0f;
    pEM = micros();

    for (int i = 0; i < Config::VHS; i++) vH[i] = {0.0f, 0.0f, 0UL};
    for (int i = 0; i < Config::AHS; i++) aH[i] = {0.0f, 0.0f, 0UL};
    for (int i = 0; i < Config::PHS; i++) pH[i] = {0.0f, 0.0f, 0.0f, 0UL};

    Serial.println("Enhanced system ready. Rotate motor to begin estimation.");
}

void loop() {
    eSC();
    float cRAK = gIEA();
    uUEA(cRAK);
    aEKF(cRAK);

    eMS = round(uEA / Config::TPF * Config::ESR);
    
    long newRC = round(uEA / Config::TPF);
    if (abs(newRC - rC) < 1000000) { // Overflow protection
        rC = newRC;
    }

    lD();

    if (Serial.available() > 0) {
        char cmd = toupper(Serial.read());
        while (Serial.available() > 0) Serial.read();
        if (cmd == 'T') {
            Serial.println("\nRe-tuning requested...");
            aTEP();
            
            int iNS = Config::DNS;
            long rI[3] = {0};
            for (int i = 0; i < iNS; i++) {
                pinMode(Config::CP, OUTPUT);
                digitalWrite(Config::CP, LOW);
                TCNT1 = 0;
                analogWrite(Config::CP, EPW);
                delayMicroseconds(Config::EPWU);
                analogWrite(Config::CP, 0);
                pinMode(Config::CP, INPUT);
                delayMicroseconds(Config::FSTU);
                rI[0] += cAR(Config::PA);
                rI[1] += cAR(Config::PB);
                rI[2] += cAR(Config::PC);
            }
            
            // BUG FIX: Perform floating point division to avoid precision loss.
            BLA::Matrix<NKF1, 1> avg_r;
            for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(rI[i]) / iNS;
            aMKF(avg_r);

            float nIRA = gIEA();
            xe(0) = nIRA;
            xe(1) = 0.0f;
            Pe.Fill(0.0f);
            Pe(0, 0) = 1.0f;
            Pe(1, 1) = 1.0f;
            pEM = micros();

            fAR = true;
            uUEA(nIRA);
            
            rC = 0;
            eMS = 0;

            vHI = aHI = pHI = 0;
            vHF = aHF = pHF = false;
            
            sV = vT = 0.0f;
            cDC = 0;
            iS = false;
            iCP = false;
            eL = false;

            Serial.println("Re-tuning complete. Enhanced system ready.");
        }
    }
}

void initM() {
    Ik1.Fill(0.0f);
    Ie.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Ik1(i, i) = 1.0f;
    for (int i = 0; i < NEKF; i++) Ie(i, i) = 1.0f;

    Qk1.Fill(0.0f);
    Rk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) {
        Qk1(i, i) = Config::QDK1;
        Rk1(i, i) = Config::RDK1;
    }

    Qe.Fill(0.0f);
    Qe(0, 0) = Config::QAE;
    Qe(1, 1) = Config::QVE;
    Re(0, 0) = Config::RAE;
}

void aTEP() {
    Serial.println("\n--- Starting Auto-Tune ---");
    Serial.println("Ensure motor is stationary!");

    xk1.Fill(0.0f);
    Pk1.Fill(0.0f);
    for (int i = 0; i < NKF1; i++) Pk1(i, i) = 1000.0f;

    for (int iter = 0; iter < Config::PCI; iter++) {
        int aTS = Config::DNS;
        long rT[3] = {0};
        for (int i = 0; i < aTS; i++) {
            pinMode(Config::CP, OUTPUT);
            digitalWrite(Config::CP, LOW);
            TCNT1 = 0;
            analogWrite(Config::CP, EPW);
            delayMicroseconds(Config::EPWU);
            analogWrite(Config::CP, 0);
            pinMode(Config::CP, INPUT);
            delayMicroseconds(Config::FSTU);
            rT[0] += cAR(Config::PA);
            rT[1] += cAR(Config::PB);
            rT[2] += cAR(Config::PC);
        }
        
        // BUG FIX: Perform floating point division to avoid precision loss.
        BLA::Matrix<NKF1, 1> avg_rT;
        for (int i = 0; i < 3; i++) avg_rT(i) = static_cast<float>(rT[i]) / aTS;
        aMKF(avg_rT);
        
        int cAA = static_cast<int>((xk1(0) + xk1(1) + xk1(2)) / 3.0f);
        int err = cAA - Config::TQA;

        if (abs(err) <= Config::PCT) {
            Serial.println("--- Auto-tune complete! ---");
            break;
        }

        EPW += (err < 0) ? Config::PAS : -Config::PAS;
        EPW = constrain(EPW, 1, 255);

        delay(5);
    }
    Serial.print("Final Tuned PWM Value: "); Serial.println(EPW);
}

int gDNS(float cV) {
    float aV = fabsf(cV);
    
    if (aV < Config::LV_TH) {
        return Config::MAX_S;
    } else if (aV > Config::HV_TH) {
        return Config::MIN_S;
    } else {
        float nV = (aV - Config::LV_TH) / (Config::HV_TH - Config::LV_TH);
        int nS = Config::MAX_S - static_cast<int>(nV * (Config::MAX_S - Config::MIN_S));
        return constrain(nS, Config::MIN_S, Config::MAX_S);
    }
}

void eSC() {
    long r[3] = {0};
    
    float cV = fabsf(xe(1));
    int nST = gDNS(cV);

    for (int i = 0; i < nST; i++) {
        pinMode(Config::CP, OUTPUT);
        digitalWrite(Config::CP, LOW);
        TCNT1 = 0;
        analogWrite(Config::CP, EPW);
        delayMicroseconds(Config::EPWU);
        analogWrite(Config::CP, 0);
        pinMode(Config::CP, INPUT);
        delayMicroseconds(Config::FSTU);
        r[0] += cAR(Config::PA);
        r[1] += cAR(Config::PB);
        r[2] += cAR(Config::PC);
    }
    
    // BUG FIX: Perform floating point division to avoid precision loss.
    BLA::Matrix<NKF1, 1> avg_r;
    for (int i = 0; i < 3; i++) avg_r(i) = static_cast<float>(r[i]) / nST;
    aMKF(avg_r);
    
    pinMode(Config::CP, OUTPUT);
    digitalWrite(Config::CP, LOW);
}

// BUG FIX: Corrected the Kalman Filter algorithm implementation.
void aMKF(const BLA::Matrix<NKF1, 1>& z) { // readings
    // --- Predict Step ---
    // State prediction is xk1 = xk1 (no change model)
    // Covariance prediction
    Pk1 += Qk1;

    // --- Update Step ---
    BLA::Matrix<NKF1, 1> y = z - xk1; // Innovation
    BLA::Matrix<NKF1, NKF1> S = Pk1 + Rk1; // Innovation covariance

    // Numerical stability check
    float det = BLA::Determinant(S);
    if (fabsf(det) < 1e-10f) {
        Serial.println("Warning: KF Matrix S near singular, skipping update");
        return;
    }
    
    auto K = Pk1 * BLA::Inverse(S); // Kalman Gain
    xk1 += K * y; // Update state estimate
    Pk1 = (Ik1 - K) * Pk1; // Update covariance estimate

    uPH(xk1(0), xk1(1), xk1(2), micros());
}

float gIEA() {
    float a = xk1(0), b = xk1(1), c = xk1(2);
    a -= Config::TQA;
    b -= Config::TQA;
    c -= Config::TQA;
    
    float al = (2.0f / 3.0f) * (a - 0.5f * (b + c)); // alpha
    float be = (1.0f / Config::SQ3) * (b - c); // beta (Original code had 2/3 factor, standard Clarke is 1/sqrt(3))
    return atan2f(be, al);
}

void uUEA(float cREA) {
    if (fAR) {
        uEA = cREA;
        pREA = cREA;
        fAR = false;
        return;
    }

    float aD = cREA - pREA;

    if (aD > Config::PF) aD -= Config::TPF;
    else if (aD < -Config::PF) aD += Config::TPF;

    uEA += aD;
    pREA = cREA;
}

void uVH(float a, float v, unsigned long ts) {
    vH[vHI] = {a, v, ts};
    vHI = (vHI + 1) % Config::VHS;
    if (vHI == 0) vHF = true;
}

void uAH(float a, float v, unsigned long ts) {
    aH[aHI] = {a, v, ts};
    aHI = (aHI + 1) % Config::AHS;
    if (aHI == 0) aHF = true;
}

void uPH(float pA, float pB, float pC, unsigned long ts) {
    pH[pHI] = {pA, pB, pC, ts};
    pHI = (pHI + 1) % Config::PHS;
    if (pHI == 0) pHF = true;
}

float gSV() {
    if (!vHF && vHI < 3) return xe(1);
    
    float sK = 0.0f, cK = 0.0f;
    float wS = 0.0f;
    int count = vHF ? Config::VHS : vHI;
    
    for (int i = 0; i < count; i++) {
        int idx = (vHI - 1 - i + Config::VHS) % Config::VHS;
        float w = exp(-i * 0.2f);
        kS(vH[idx].velocity * w, sK, cK);
        wS += w;
    }
    
    return (wS > 0) ? sK / wS : xe(1);
}

float gVT() {
    if (!vHF && vHI < 5) return 0.0f;
    
    int count = min(vHF ? Config::VHS : vHI, 5);
    float sXk = 0.0f, cX = 0.0f;
    float sYk = 0.0f, cY = 0.0f;
    float sXYk = 0.0f, cXY = 0.0f;
    float sX2k = 0.0f, cX2 = 0.0f;
    
    for (int i = 0; i < count; i++) {
        int idx = (vHI - 1 - i + Config::VHS) % Config::VHS;
        float x = i;
        float y = vH[idx].velocity;
        kS(x, sXk, cX);
        kS(y, sYk, cY);
        kS(x * y, sXYk, cXY);
        kS(x * x, sX2k, cX2);
    }
    
    float den = count * sX2k - sXk * sXk;
    return (fabsf(den) > 1e-6f) ? (count * sXYk - sXk * sYk) / den : 0.0f;
}

bool iVC() {
    if (!vHF && vHI < 3) return false;
    
    int count = min(vHF ? Config::VHS : vHI, 5);
    
    float meanK = 0.0f, cM = 0.0f;
    for (int i = 0; i < count; i++) {
        int idx = (vHI - 1 - i + Config::VHS) % Config::VHS;
        kS(vH[idx].velocity, meanK, cM);
    }
    float mean = meanK / count;
    
    float vK = 0.0f, cV = 0.0f;
    
    for (int i = 0; i < count; i++) {
        int idx = (vHI - 1 - i + Config::VHS) % Config::VHS;
        float diff = vH[idx].velocity - mean;
        kS(diff * diff, vK, cV);
    }
    vK /= count;
    
    return sqrt(vK) < Config::VNT;
}

bool cPC() {
    if (!pHF && pHI < 5) return false;
    
    int count = pHF ? Config::PHS : pHI;
    
    float sAk = 0.0f, cA = 0.0f;
    
    count = min(count, Config::PHS);
    float a[Config::PHS];

    for (int i = 0; i < count; ++i) {
        int idx = (pHI - 1 - i + Config::PHS) % Config::PHS;
        float aP = pH[idx].phaseA;
        float bP = pH[idx].phaseB;
        float cP = pH[idx].phaseC;
        
        aP -= Config::TQA;
        bP -= Config::TQA;
        cP -= Config::TQA;

        float al = (2.0f / 3.0f) * (aP - 0.5f * (bP + cP));
        float be = (1.0f / Config::SQ3) * (bP - cP); // Corrected Clarke transform
        a[i] = atan2f(be, al);
        kS(a[i], sAk, cA);
    }

    float mA = sAk / count;
    float vK = 0.0f, cV = 0.0f;
    for (int i = 0; i < count; ++i) {
        float diff = nA(a[i] - mA);
        kS(diff * diff, vK, cV);
    }
    vK /= count;
    
    return sqrt(vK) < Config::AHNT;
}

bool iRS() {
    float trend = gVT();
    float cV = fabsf(sV);
    
    return (cV < Config::MVC && trend < -0.01f && iVC());
}

void aEKF(float mA) {
    unsigned long now = micros();
    float dt = (now - pEM) / 1e6f;
    if (dt <= 0) dt = 1e-6; // Prevent division by zero or negative time
    pEM = now;

    float angle = xe(0), velocity = xe(1);

    float pA = nA(angle + velocity * dt);
    float pV = velocity;

    BLA::Matrix<NEKF, 1> xP;
    xP(0) = pA;
    xP(1) = pV;

    BLA::Matrix<NEKF, NEKF> A;
    A.Fill(0.0f);
    A(0, 0) = 1.0f; A(0, 1) = dt;
    A(1, 1) = 1.0f;

    Pe = A * Pe * A.Transpose() + Qe;

    float y = nA(mA - pA);
    BLA::Matrix<1, 1> S = Re + Pe(0, 0);
    
    if (fabsf(S(0, 0)) < 1e-10f) {
        Serial.println("Warning: EKF S matrix near singular");
        return;
    }
    
    BLA::Matrix<NEKF, 1> K;
    K(0) = Pe(0, 0) / S(0, 0);
    K(1) = Pe(1, 0) / S(0, 0);
    
    xe = xP + K * y;
    xe(0) = nA(xe(0)); // Normalize angle after update

    BLA::Matrix<NEKF, NEKF> IKH = Ie;
    IKH(0, 0) -= K(0);
    IKH(1, 0) -= K(1);
    Pe = IKH * Pe;

    uVH(xe(0), xe(1), now);
    uAH(xe(0), xe(1), now);
    
    sV = gSV();
    vT = gVT();
    iS = iRS();
    iCP = cPC();
    bool iCV = iVC();
    float vM = fabsf(sV);
    bool iM = vM > Config::VST;

    if (iM && iCV && iCP && !iS) {
        if (!eL && cDC > 3) {
            eL = true;
            Serial.println("EKF locked - confident movement detected.");
        }
        lNZM = now;
    } else if (eL) {
        bool sU = false;
        if (iS) {
            sU = true;
            Serial.println("EKF unlocked - rotation stopping detected.");
        } else if (!iM && (now - lNZM > Config::ZVLT)) {
            sU = true;
            Serial.println("EKF unlocked - prolonged zero velocity.");
        } else if (!iCV && vM < Config::MVC) {
            sU = true;
            Serial.println("EKF unlocked - inconsistent low velocity.");
        } else if (!iCP) {
            sU = true;
            Serial.println("EKF unlocked - inconsistent phase signals.");
        }
        if (sU) {
            eL = false;
            cDC = 0;
        }
    }
}

float nA(float a) {
    while (a <= -Config::PF) a += Config::TPF;
    while (a > Config::PF) a -= Config::TPF;
    return a;
}

int cAR(int p) {
    int val = analogRead(p);
    // analogRead on AVR is guaranteed to be 0-1023, so checks are redundant but harmless.
    return val; 
}

void lD() {
    Serial.print(xk1(0)-Config::TQA); Serial.print(",");
    Serial.print(xk1(1)-Config::TQA); Serial.print(",");
    Serial.print(xk1(2)-Config::TQA); Serial.print(",");
    Serial.print(uEA, 3); Serial.print(",");
    Serial.print(sV, 4); Serial.print(",");
    Serial.print(eL ? 1 : 0); Serial.print(",");
    Serial.print(iCP ? 1 : 0);Serial.print(",");
    Serial.print(eMS);Serial.print(",");
    Serial.println(rC);
} 
