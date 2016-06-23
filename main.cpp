#include <iostream>
#include <ctime>

#include "Simulator.h"

int main()

{
    bool writeLog = false;

    std::clock_t begin = std::clock();

    // Sizes are in micrometers (μm)
    double d = 4;
	double h = 1;

    // Time in seconds. For different d, use different t_s.
    // d = 4μm : Symbol duration (t_s) = 0.213 s
    // d = 8μm : Symbol duration (t_s) = 0.949 s
    // d = 16μm : Symbol duration (t_s) = 4.064 s
    // d = 32μm : Symbol duration (t_s) = 17.391 s
    double t_s = 0.213;

	double r_cell = 5;
	double r_molecule = 0.0025;

	double D = 79.4;
	double TIME_STEP = 0.001;
	double SIGMA = sqrt(2 * D * TIME_STEP);
	double NUM_MOLS = 100000;
	
    Simulator mysimulator(d, h, t_s, r_cell, r_molecule, NUM_MOLS, SIGMA, TIME_STEP, writeLog);
	mysimulator.runSimulation();


    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << elapsed_secs << std::endl;

	return 0;
}
