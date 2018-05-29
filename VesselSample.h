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
    
    void compute(list<Harmonic> & in);
private:
    const static double CONST_VISCOSITY = 5000000.0;
    const static double CONST_DENSITY   = 1060.0;

    const static double YOUNG_K1 = 8000000.0;
    const static double YOUNG_K2 = -2253.0;
    const static double YOUNG_K3 = 865000.0;

    const static double LAPLACE_A1 = 4.0;
    const static double LAPLACE_A2 = 50.0;
    const static double LAPLACE_B1 = 6.0;
    const static double LAPLACE_B2 = 70.0;
    
    double radius;
    double length;
    double wall;
    
    complex<double> getF10(double womersly);
    double getWomersley(double frequency);
    double getYoungModulus();
    double getStaticYoungModulus(double frequency);
    double getR1(double frequency);
    double getR2(double frequency);
    double getR3(double frequency);
    double getL1(double frequency);
    double getC1(double frequency);
    double getC2(double frequency);
    double getC3(double frequency);
    complex<double> getZx(double frequency);
    complex<double> getZy(double frequency);
    complex<double> getZforParallelRC(double r, double c, double w);
    void processHarmonic(Harmonic & harmonic);
};

#endif /* VESSELSAMPLE_H */

