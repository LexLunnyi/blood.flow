#ifndef VESSELSAMPLE_H
#define VESSELSAMPLE_H

#include "Harmonic.h"

#include <cmath>
#include <complex>
#include <math.h>
#include <list>

using namespace std;


class VesselSample {
public:
    VesselSample(double radius, double length, double wall);
    VesselSample(const VesselSample& orig);
    virtual ~VesselSample();
    
    void compute(list<PHarmonic> & in);
    void processHarmonic(PHarmonic & harmonic);
    void setFrequency(double frequency);
private:
    bool wasDebug = false;
    
    const double CONST_VISCOSITY = 5000000.0;
    const double CONST_DENSITY   = 1.060;

    const double YOUNG_K1 = 8000000.0;
    const double YOUNG_K2 = -2253.0;
    const double YOUNG_K3 = 865000.0;

    const double LAPLACE_A1 = 4.0;
    const double LAPLACE_A2 = 50.0;
    const double LAPLACE_B1 = 6.0;
    const double LAPLACE_B2 = 70.0;
    
    double radius;
    double length;
    double wall;
    double young;
    
    double angularVelocity;
    double womersleyNumber;
    complex <double> f10;
    double staticYoung;
    
    double R1;
    double L1;
    double C1;
    double R2;
    double C2;
    double R3;
    double C3;
    
    complex<double> Zx;
    complex<double> Zy;
    
    complex<double> getF10();
    double getWomersley();
    double getYoungModulus();
    double getStaticYoungModulus();
    double getR1();
    double getR2();
    double getR3();
    double getL1();
    double getC1();
    double getC2();
    double getC3();
    complex<double> getZx();
    complex<double> getZy();
    complex<double> getZforParallelRC(double r, double c, double w);
    
    complex<double> complexBessel(unsigned int order, complex<double> in);
};

#endif /* VESSELSAMPLE_H */