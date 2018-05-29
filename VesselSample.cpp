#include <algorithm>

#include "VesselSample.h"



VesselSample::VesselSample(double radius, double length, double wall) {
    this->radius = radius;
    this->length = length;
    this->wall = wall;
}

VesselSample::VesselSample(const VesselSample& orig) {
}

VesselSample::~VesselSample() {
}


double VesselSample::getWomersley(double frequency) {
    double w = 2*M_PI*frequency;
    double v = CONST_DENSITY / CONST_VISCOSITY;
    return radius*sqrt(w/v);
}



complex<double> VesselSample::getF10(double womersly) {
    complex<double> j(0, 1);
    complex<double> two(2, 0);
    complex<double> besArg = pow(j, 1.5)*womersly;
    
    complex<double> up (j1(besArg.real()), j1(besArg.imag()));
    complex<double> down (j0(besArg.real()), j0(besArg.imag()));

    return (two*up)/(besArg*down);
}


/*
 METTE S. OLUFSEN
 Structured tree outflow condition for blood flow in larger systemic arteries, 1997
 Formula (6) Figure 3
 r = 0 - 0.8 cm
 */
double VesselSample::getYoungModulus() {
    return (YOUNG_K1*exp(radius*YOUNG_K2) + YOUNG_K3);
}



double VesselSample::getStaticYoungModulus(double frequency) {
    double w = 2*M_PI*frequency;
    complex<double> a1w(LAPLACE_A1, w);
    complex<double> a2w(LAPLACE_A2, w);
    complex<double> b1w(LAPLACE_B1, w);
    complex<double> b2w(LAPLACE_B2, w);
    complex<double> a1(LAPLACE_A1, 0);
    complex<double> a2(LAPLACE_A2, 0);
    complex<double> b1(LAPLACE_B1, 0);
    complex<double> b2(LAPLACE_B2, 0);
    complex<double> young(getYoungModulus(), 0);
    complex<double> koeff = (b1*b2*a1w*a2w)/(a1*a2*b1w*b2w);
    return (young / koeff).real();
}



double VesselSample::getR1(double frequency) {
    double w = 2*M_PI*frequency;
    
    complex<double> one(1, 0);
    complex<double> koeff(CONST_DENSITY/(M_PI*pow(radius, 2)), 0);
    complex<double> cw(-1*w, 0);
        
    complex<double> r1 = koeff*(cw/(one - getF10(getWomersley(frequency))));
    return r1.imag();
}



double VesselSample::getL1(double frequency) {
    double w = 2*M_PI*frequency;
    
    complex<double> one(1, 0);
    complex<double> koeff(CONST_DENSITY/(M_PI*radius*radius), 0);
        
    complex<double> l1 = koeff*(one/(one - getF10(getWomersley(frequency))));
    return l1.real();
}



double VesselSample::getC1(double frequency) {
    double up = 3*M_PI*pow(radius, 2)*pow(radius + wall, 2);
    double down = wall*(2*radius + wall)*getStaticYoungModulus(frequency);
    return up / down;
}




double VesselSample::getC2(double frequency) {
    double up = LAPLACE_A1*LAPLACE_A2*(LAPLACE_B2-LAPLACE_B1);
    double down = LAPLACE_B2*(LAPLACE_B1-LAPLACE_A1)*(LAPLACE_A2-LAPLACE_B1);
    return (up/down)*getC1(frequency);
}



double VesselSample::getC3(double frequency) {
    double up = LAPLACE_A1*LAPLACE_A2*(LAPLACE_B2-LAPLACE_B1);
    double down = LAPLACE_B1*(LAPLACE_B2-LAPLACE_A1)*(LAPLACE_B2-LAPLACE_A2);
    return (up/down)*getC1(frequency);
}




double VesselSample::getR2(double frequency) {
    return 1.0/(LAPLACE_B1*getC2(frequency));
}



double VesselSample::getR3(double frequency) {
    return 1.0/(LAPLACE_B2*getC3(frequency));
}




complex<double> VesselSample::getZx(double frequency) {
    double w = 2*M_PI*frequency;
    complex<double> res(getR1(frequency), w*getL1(frequency));
    return res;
}



complex<double> VesselSample::getZy(double frequency) {
    complex<double> res(0, 0);
    
    double w = 2*M_PI*frequency;
    complex<double> z1(0, -1.0/(w*getC1(frequency)));
    complex<double> z2 = getZforParallelRC(getR2(frequency), getC2(frequency), w);
    complex<double> z3 = getZforParallelRC(getR3(frequency), getC3(frequency), w);
    
    return z1 + z2 + z3;
}



complex<double> VesselSample::getZforParallelRC(double r, double c, double w) {
    complex<double> res(0, 0);
    
    double up = r*pow(w,2)*pow(c,2);
    double down = pow(r,2) + pow(w,2)*pow(c,2);

    res.real(up/down);
    
    up = pow(r,2)*w*c;
    down = pow(r,2) + pow(w,2)*pow(c,2);
    
    res.imag(up/down);
    
    return res;
}




void VesselSample::compute(list<Harmonic> & in) {
    for (auto it = in.begin(); it != in.end(); ++it) {
        processHarmonic(*it);
    }
}



void VesselSample::processHarmonic(Harmonic & harmonic) {
    complex<double> Uin(harmonic.inAmplitude*cos(harmonic.inPhase), harmonic.inAmplitude*sin(harmonic.inPhase));
    complex<double> zx = getZx(harmonic.angularVelocity);
    complex<double> zy = getZy(harmonic.angularVelocity);
    complex<double> Uout = (Uin*zy)/(zx+zy);
    
    harmonic.outAmplitude = sqrt(pow(Uout.real(), 2) + pow(Uout.imag(),2));
    harmonic.outPhase = atan(Uout.real()/Uout.imag());
}