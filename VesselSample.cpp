#include <algorithm>

#include "VesselSample.h"



VesselSample::VesselSample(double radius, double length, double wall) {
    this->radius = radius;
    this->length = length;
    this->wall = wall;
    this->young = getYoungModulus();
}

VesselSample::VesselSample(const VesselSample& orig) {
}

VesselSample::~VesselSample() {
}



void VesselSample::setFrequency(double frequency) {
    //this->angularVelocity = 2.0*M_PI*frequency;
    this->angularVelocity = frequency/(2.0*M_PI);
    this->womersleyNumber = getWomersley();
    this->f10 = getF10();
    this->staticYoung = getStaticYoungModulus();
    this->R1 = getR1();
    this->L1 = getL1();
    this->C1 = getC1();
    this->C2 = getC2();
    this->C3 = getC3();
    this->R2 = getR2();
    this->R3 = getR3();
    this->Zx = getZx();
    this->Zy = getZy();
}


double VesselSample::getWomersley() {
    double v = CONST_DENSITY / CONST_VISCOSITY;
    return radius*sqrt(this->angularVelocity/v);
}

 
/*
complex<double> VesselSample::getF10() {
    complex<double> j(0, 1);
    complex<double> two(2, 0);
    complex<double> complexWomersley(this->womersleyNumber, 0);
    complex<double> besArg = pow(j, 1.5)*complexWomersley;
    
    complex<double> up = complexBessel(1, besArg);
    complex<double> down = complexBessel(0, besArg);
    
    //complex<double> up (j1(besArg.real()), j1(besArg.imag()));
    //complex<double> down (j0(besArg.real()), j0(besArg.imag()));

    return (two*up)/(besArg*down);
}
 */


/*
 METTE S. OLUFSEN
 Structured tree outflow condition for blood flow in larger systemic arteries, 1997
 Formula (6) Figure 3
 r = 0 - 0.8 cm
 */
double VesselSample::getYoungModulus() {
    return (YOUNG_K1*exp(radius*YOUNG_K2) + YOUNG_K3);
}



double VesselSample::getStaticYoungModulus() {
    complex<double> a1w(LAPLACE_A1, this->angularVelocity);
    complex<double> a2w(LAPLACE_A2, this->angularVelocity);
    complex<double> b1w(LAPLACE_B1, this->angularVelocity);
    complex<double> b2w(LAPLACE_B2, this->angularVelocity);
    complex<double> a1(LAPLACE_A1, 0);
    complex<double> a2(LAPLACE_A2, 0);
    complex<double> b1(LAPLACE_B1, 0);
    complex<double> b2(LAPLACE_B2, 0);
    complex<double> young(this->young, 0);
    complex<double> koeff = (b1*b2*a1w*a2w)/(a1*a2*b1w*b2w);
    return (young / koeff).real();
}



double VesselSample::getR1() {   
    complex<double> one(1, 0);
    complex<double> koeff(CONST_DENSITY/(M_PI*pow(radius, 2)), 0);
    complex<double> cw(-1*this->angularVelocity, 0);
        
    complex<double> r1 = koeff*(cw/(one - this->f10));
    return r1.imag();
}



double VesselSample::getL1() {
    complex<double> one(1, 0);
    complex<double> koeff(CONST_DENSITY/(M_PI*pow(radius, 2)), 0);
    complex<double> cw(this->angularVelocity, 0);
        
    complex<double> l1 = koeff*(one/(one - this->f10));
    return l1.real();
}



double VesselSample::getC1() {
    double up = 3*M_PI*pow(radius, 2)*pow(radius + wall, 2);
    double down = wall*(2*radius + wall)*this->staticYoung;
    return up / down;
}




double VesselSample::getC2() {
    double up = LAPLACE_A1*LAPLACE_A2*(LAPLACE_B2-LAPLACE_B1);
    double down = LAPLACE_B2*(LAPLACE_B1-LAPLACE_A1)*(LAPLACE_A2-LAPLACE_B1);
    return (up/down)*this->C1;
}



double VesselSample::getC3() {
    double up = LAPLACE_A1*LAPLACE_A2*(LAPLACE_B2-LAPLACE_B1);
    double down = LAPLACE_B1*(LAPLACE_B2-LAPLACE_A1)*(LAPLACE_B2-LAPLACE_A2);
    return (up/down)*this->C1;
}




double VesselSample::getR2() {
    return 1.0/(LAPLACE_B1*this->C2);
}



double VesselSample::getR3() {
    return 1.0/(LAPLACE_B2*this->C3);
}




complex<double> VesselSample::getZx() {
    complex<double> res(this->R1, this->angularVelocity*this->L1);
    return res;
}



complex<double> VesselSample::getZy() {
    complex<double> z1(0, -1.0/(this->angularVelocity*this->C1));
    complex<double> z2 = getZforParallelRC(this->R2, this->C2, this->angularVelocity);
    complex<double> z3 = getZforParallelRC(this->R3, this->C3, this->angularVelocity);
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




void VesselSample::compute(list<PHarmonic> & in) {
    list<PHarmonic>::iterator it;
    for (it = in.begin(); it != in.end(); ++it) {
        processHarmonic(*it);
        if (!wasDebug) {
            cout << "FREQUENCY: " << (*it)->frequency <<  "; Womersley: " << womersleyNumber << "; F10: " << f10 << "; R1: " << R1 << "; L1: " << L1 << "; C1: " << C1 << "; R2: " << R2 << "; C2: " << C2 << "; R3: " << R3 << "; C3: " << C3 << endl;   
        }
    }
    if (!wasDebug) {
        wasDebug = true;
    }
}



void VesselSample::processHarmonic(PHarmonic & harmonic) {
    setFrequency(harmonic->frequency);
    complex<double> Uin(harmonic->inAmplitude*cos(harmonic->inPhase), harmonic->inAmplitude*sin(harmonic->inPhase));
    complex<double> Uout = (Uin*this->Zy)/(this->Zx + this->Zy);
    
    harmonic->outAmplitude = sqrt(pow(Uout.real(), 2) + pow(Uout.imag(),2));
    harmonic->outPhase = atan(Uout.real()/Uout.imag());
}



complex<double> VesselSample::complexBessel(unsigned int order, complex<double> in) {
    complex<double> res;
    if (0 == order) {
        res.real(j0(in.real())+y0(in.imag()));
        res.imag(0);
    } else {
        res.real(j1(in.real()));
        res.imag(y1(in.imag()));        
    }
    return res;
}


complex<double> VesselSample::getF10() {
    complex<double> betta(0, (pow(this->womersleyNumber, 2)/16));
    complex<double> one(1, 0);
    complex<double> two(2, 0);
    
    complex<double> up = sqrt(one + betta);
    complex<double> down = sqrt(one + betta) + two*betta;
    
    return up/down;
}