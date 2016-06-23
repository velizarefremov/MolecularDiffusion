#ifndef INC_SIMULATOR_H
#define INC_SIMULATOR_H

#include <list>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

#include "Transmitter.h"
#include "Receiver.h"
#include "Molecule.h"


class Simulator
{
public:

    Simulator(double _d, double _h, double _t_s, double _r_cell, double _r_molecule, int _numMolecules, double _SIGMA, double _TIME_STEP, bool _writeLog)
	{
		numMolecules = _numMolecules;
		SIGMA = _SIGMA;
		TIME_STEP = _TIME_STEP;

		currentTime = 0;
		d = _d;
		h = _h;
		t_s = _t_s;
		r_cell = _r_cell;
		r_molecule = _r_molecule;

		r_environment = 100 * r_cell;

		r_envcoll = r_environment - r_molecule;

		lostMolCount = 0;

        writeLog = _writeLog;

        if(writeLog)
        {
            logstring.flush();
            // logstring << new SimpleDateFormat("dd-MM-yyyy h:mm a").format(new Date());
            logstring << std::endl;
        }

        wrongprob[0] = 0;
        wrongprob[1] = 0;
        wrongprob[2] = 0;
        wrongprob[3] = 0;

		initDevices();

		initMolecules();
	}

	void initMolecules()
	{
		molecules.clear();

		for (std::list<Transmitter>::iterator it = transmitters.begin(); it != transmitters.end(); it++)
		{
			// cout << "HERE" << endl;

			Transmitter tr = *it;

			int destin = tr.getDestination();
			double xpos = tr.getX() + tr.getRadius();
			double ypos = tr.getY();
			double zpos = tr.getZ();

            std::vector<Molecule> toAdd(numMolecules, Molecule(xpos, ypos, zpos, destin));

			molecules.insert(molecules.end(), toAdd.begin(), toAdd.end());

			// cout << toAdd.size() << endl;
		}

        std::cout << molecules.size() << std::endl;
	}

	void initDevices()
	{
		// Change this part for device positions.

		// Assume R2 is in the center for now.

		double single_dist = 2 * r_cell + h;
		double single_trans_dist = 2 * r_cell + d;

		receivers.push_front(Receiver(0, -2 * single_dist, 0, r_cell, 4, r_molecule)); // R4
		receivers.push_front(Receiver(0, -single_dist, 0, r_cell, 3, r_molecule)); // R3
		receivers.push_front(Receiver(0, 0, 0, r_cell, 2, r_molecule)); // R2
		receivers.push_front(Receiver(0, single_dist, 0, r_cell, 1, r_molecule)); // R1

		transmitters.push_front(Transmitter(-single_trans_dist, -2 * single_dist, 0, r_cell, 4, 4, r_molecule)); // T4
		transmitters.push_front(Transmitter(-single_trans_dist, -single_dist, 0, r_cell, 3, 3, r_molecule)); // T3
		transmitters.push_front(Transmitter(-single_trans_dist, 0, 0, r_cell, 2, 2, r_molecule)); // T2
		transmitters.push_front(Transmitter(-single_trans_dist, single_dist, 0, r_cell, 1, 1, r_molecule)); // T1

        std::cout << receivers.size() << std::endl;
        std::cout << transmitters.size() << std::endl;

        if(writeLog)
        {
            logstring << "T1;" << this->d << ";" << this->h << ";" << this->r_cell << ";" << this->r_molecule << ";" << 79.4 << ";" << this->numMolecules << std::endl;
            logstring << "T2;" << this->d << ";" << this->h << ";" << this->r_cell << ";" << this->r_molecule << ";" << 79.4 << ";" << this->numMolecules << std::endl;
            logstring << "T3;" << this->d << ";" << this->h << ";" << this->r_cell << ";" << this->r_molecule << ";" << 79.4 << ";" << this->numMolecules << std::endl;
            logstring << "T4;" << this->d << ";" << this->h << ";" << this->r_cell << ";" << this->r_molecule << ";" << 79.4 << ";" << this->numMolecules << std::endl;
        }
	}

    double getDistance(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
    }

    double gaussian(double mu, double sigma)
    {
        return mu + sigma * distribution(generator);
    }

	void advanceTime() {
		
		int totreceivers = receivers.size();

		for (std::vector<Molecule>::iterator it = molecules.begin(); it != molecules.end(); ++it)
		{
			if (it->getDestination() != -1) {

                updateLocation(transmitters, SIGMA, *it);

				// Check if molecule went out of the simulation radius which is 100*r_cell

				if (it->checkBoundary(r_envcoll))
				{
					it->setDestination(-1);
					lostMolCount++;
				}
				else
				{
					for (std::list<Receiver>::iterator it2 = receivers.begin(); it2 != receivers.end(); ++it2)
					{
                        if (isMoleculeReceived(*it, *it2)) {
							it->setDestination(-1);
							break;
						}
					}
				}
			}
		}

        // Log File
        // logstring << df.format(currentTime + this->TIME_STEP) + ";";

        int receiver_index = 0;
        for (std::list<Receiver>::iterator it3 = receivers.begin(); it3 != receivers.end(); ++it3)
        {
            logstring << it3->getNumberOfMoleculesReceived();

            if ((receiver_index + 1) != totreceivers)
            {
                logstring << ";";
			}

            receiver_index++;
		}

        logstring << std::endl;

	}

	void runSimulation()
	{
		// Simulate up to (2 * t_s)

		while (currentTime <= t_s) {
			advanceTime();
			currentTime += TIME_STEP;

		}

		// Here get result at t_s
		print();


		// Now to 2 * t_s
		while (currentTime <= 2 * t_s) {
			advanceTime();
			currentTime += TIME_STEP;
		}

		// Here get result at 2*t_s    
		print();

        if(writeLog)
        {
            printfile();
        }
	}

	void print()
	{
        std::cout << "Number of total molecules sent foreach receiver: " << numMolecules << std::endl;
        std::cout << "Number of total molecules lost at the boundary: " << lostMolCount << std::endl;

		int totreceivers = receivers.size();

		for (std::list<Receiver>::iterator it2 = receivers.begin(); it2 != receivers.end(); ++it2)
		{
            std::cout << "Number of molecules received R" + it2->getID() << ": " << it2->getNumberOfMoleculesReceived() << std::endl;
            std::cout << "Number of molecules wrongly received R" << it2->getID() << ": " << it2->getNumberOfFalseMoleculesReceived() << std::endl;
			
            std::cout << "Probability to capture (simulation) R" << it2->getID() << ": " << it2->getNumberOfMoleculesReceived() / (double)numMolecules << std::endl;
            std::cout << "Probability to wrong capture (simulation) R" << it2->getID() << ": " << it2->getNumberOfFalseMoleculesReceived() / (double)numMolecules << std::endl;
            std::cout << std::endl;
		}

        std::cout << std::endl;

        std::cout << wrongprob[0] << std::endl;
        std::cout << wrongprob[1] << std::endl;
        std::cout << wrongprob[2] << std::endl;
        std::cout << wrongprob[3] << std::endl;
	}

    // Molecule Location
    void updateLocation(std::list<Transmitter>& transmitters, double sigma, Molecule& mol)
    {
        double tmp_x = 0;
        double tmp_y = 0;
        double tmp_z = 0;

        bool isCollided = true;

        while (isCollided)
        {
            isCollided = false;

            tmp_x = mol.getX() + gaussian(0, sigma);
            tmp_y = mol.getY() + gaussian(0, sigma);
            tmp_z = mol.getZ() + gaussian(0, sigma);

            for (std::list<Transmitter>::iterator it = transmitters.begin(); it != transmitters.end(); it++)
            {
                if (isMoleculeReceived(tmp_x, tmp_y, tmp_z, *it))
                {
                    isCollided = true;
                    break;
                }
            }
        }

        mol.setPosition(tmp_x, tmp_y, tmp_z);
    }

    // Transmitter
    bool isMoleculeReceived(double _x, double _y, double _z, Transmitter &tr)
    {
        if (getDistance(_x, _y, _z, tr.getX(), tr.getY(), tr.getZ()) <= tr.getCollradius())
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // Receiver
    bool isMoleculeReceived(Molecule& molecule, Receiver& recv)
    {
        if (getDistance(molecule.getX(), molecule.getY(), molecule.getZ(), recv.getX(), recv.getY(), recv.getZ()) <= recv.getCollradius())
        {
            int moldestid = molecule.getDestination();
            if (moldestid == recv.getID())
            {
                recv.addReceived();
            }
            else
            {
                recv.addFalseReceived();
                wrongprob[moldestid - 1]++;
            }

            return true;
        }
        else
        {
            return false;
        }
    }

	void printfile()
	{
        std::ofstream outFile;
        outFile.open("deneme.txt");
        outFile << logstring.rdbuf();
        outFile.close();
	}

private:

    std::list<Receiver> receivers;
    std::list<Transmitter> transmitters;

    std::vector<Molecule> molecules;

	double currentTime;
    int wrongprob[4];

    std::stringstream logstring;

	double d;
	double h;
	double t_s;
	double r_cell;
	double r_molecule;

	double r_environment;
	double r_envcoll;

	int lostMolCount;

	int numMolecules;
	double SIGMA;
	double TIME_STEP;

    bool writeLog;

    static std::default_random_engine generator;
    static std::normal_distribution<double> distribution;
};

std::default_random_engine Simulator::generator = std::default_random_engine(std::chrono::high_resolution_clock::now().time_since_epoch().count());
std::normal_distribution<double> Simulator::distribution = std::normal_distribution<double>(0, 1);

#endif
