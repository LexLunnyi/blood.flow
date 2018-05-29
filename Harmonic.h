#ifndef HARMONIC_H
#define HARMONIC_H

class Harmonic {
public:
    Harmonic(double amplitude, double phase, double angularVelocity);
    Harmonic(const Harmonic& orig);
    virtual ~Harmonic();
    
    double inAmplitude;
    double inPhase;
    
    double outAmplitude;
    double outPhase;
    
    double angularVelocity;
private:
};

#endif /* HARMONIC_H */

