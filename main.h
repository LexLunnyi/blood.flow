#ifndef MAIN_H
#define MAIN_H


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <memory>

#include "VesselSample.h"
#include "Harmonic.h"

using namespace std;


//#define VESSEL_RADIUS 0.0015
//#define VESSEL_LENGTH 0.0500



double inlet(double time);
double inlet_gauge(double time);

void getInputSignals(list<PHarmonic> & in);
void getGrafPoint(list<PHarmonic> & input, double time, double *pInlet, double *pOutlet);



#endif /* MAIN_H */

