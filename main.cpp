#include "main.h"

/*
 * 
 */
int main(int argc, char** argv) {
    
    //double vesselRadius = 0.005;
    //double vesselLength = 0.015;
    //double vesselWall   = 0.01;
     
    double vesselRadius = 0.0037;
    double vesselLength = 0.089;
    double vesselWall   = 0.0001;
    
    std::cout << "Test!";
    
    ofstream file("inlet.csv");
    if (!file.is_open()) {
        std::cerr << "Error open CSV file" << endl;
        return -1;
    }
    
    list<PHarmonic> input;
    getInputSignals(input);
    
    VesselSample vessel(vesselRadius, vesselLength, vesselWall);
    
    
    /*
    double frequency;
    for(int i = 1; i < 10000; i++) {
        frequency = i / 100.0;
        PHarmonic f(new Harmonic(5.0, 0.0, frequency));
        vessel.setFrequency(frequency);
        vessel.processHarmonic(f);
        file << f->outAmplitude << ";" << endl;
    }*/
    
    
    
    
    vessel.compute(input);
    double in, out = 0.0;
    for(int i = 0; i < 8000; i++) {
        getGrafPoint(input, i / 1000.0, &in, &out);
        file << 100.0+in << "; " << 80.0+out << ";" << endl;
        //file << in << "; " << out << ";" << endl;
    }
    
    
    
    file.close();
    
    
    return 0;
}




//6 + sin(x) + 5sin((sin(x)+x)) + 0.2*(1+sin(x-1.5))*sin(12.5x) +  2*(1+sin(x))*sin(x+2)
double inlet(double time) {
    double res = 6;
    res += sin(time);
    res += 5*sin(sin(time) + time);
    res += 0.2*(1+sin(time-1.5))*sin(12.5*time);
    res += 2*(1+sin(time))*sin(time+2);
    return res;
}



double inlet_gauge(double time) {
    double res = 5;
    res += 10*sin(4*time);
    res += 10*sin(6*time);
    res += sin(50*time);
    res += sin(70*time);
    return res;
}




void getInputSignals(list<PHarmonic> & in) {
    in.push_back(PHarmonic(new Harmonic(10.0, 0.0, 1.2)));
    in.push_back(PHarmonic(new Harmonic(8.0, -(M_PI/4), 2.4)));
    in.push_back(PHarmonic(new Harmonic(1.0, 0.0, 20.0)));
    //in.push_back(PHarmonic(new Harmonic(0.9, 0.0, 17.68)));
    in.push_back(PHarmonic(new Harmonic(1.0, -(M_PI/4), 50.0)));
}



void getGrafPoint(list<PHarmonic> & input, double time, double *pInlet, double *pOutlet) {
    *pInlet = 0.0;
    *pOutlet = 0.0;
    list<PHarmonic>::iterator it;
    for(it = input.begin(); it != input.end(); ++it) {
        *pInlet += (*it)->inAmplitude*sin((*it)->frequency*time + (*it)->inPhase);
        //*pOutlet += (*it)->outAmplitude*sin((*it)->frequency*time + (*it)->outPhase);
        *pOutlet += (*it)->outAmplitude*sin((*it)->frequency*time);
    }
}