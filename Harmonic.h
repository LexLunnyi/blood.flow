#ifndef HARMONIC_H
#define HARMONIC_H

#include <iostream>
#include <memory>
#include <algorithm>

using namespace std;


class Harmonic {
public:
    Harmonic(double amplitude, double phase, double frequency);
    Harmonic(const Harmonic& orig);
    virtual ~Harmonic();
    
    double inAmplitude;
    double inPhase;
    
    double outAmplitude;
    double outPhase;
    
    double frequency;
private:
};


typedef unique_ptr<Harmonic> PHarmonic;

#endif /* HARMONIC_H */