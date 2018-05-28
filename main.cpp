#include "main.h"

/*
 * 
 */
int main(int argc, char** argv) {
    double vesselRadius = 0.0025;
    
    
    std::cout << "Test!";
    
    /*
    std::complex<double> test_c(0, 1);
    std::complex<double> new_c = pow(test_c, 1.5);
    std::cout << new_c << endl;
     */
    
    ofstream file("inlet.csv");
    if (!file.is_open()) {
        std::cerr << "Error open CSV file" << endl;
        return -1;
    }
    
    double in, out = 0.0;
    for(int i = 0; i < 1000; i++) {
        in = i / 1.0;
        //out = inlet(in);
        //out = getL1(vesselRadius, in);
        file << in << "; " << out << ";" << endl;
    }
    
    file.close();
    
    
    return 0;
    
    /*
    ofstream file("inlet.csv");
    if (!file.is_open()) {
        std::cerr << "Error open CSV file" << endl;
        return -1;
    }
    
    double in, out = 0.0;
    for(int i = 0; i < 15000; i++) {
        in = i / 1000.0;
        //out = inlet(in);
        out = inlet_gauge(in);
        file << in << "; " << out << ";" << endl;
    }
    
    file.close();
    
    return 0;
     */
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