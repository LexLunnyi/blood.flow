#include "Harmonic.h"

Harmonic::Harmonic(double amplitude, double phase, double angularVelocity) {
    this->inAmplitude = amplitude;
    this->inPhase = phase;
    this->angularVelocity = angularVelocity;
    this->outAmplitude = 0.0;
    this->outPhase = 0.0;
}

Harmonic::Harmonic(const Harmonic& orig) {
}

Harmonic::~Harmonic() {
}
