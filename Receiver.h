#ifndef INC_RECEIVER_H
#define INC_RECEIVER_H

class Receiver
{
public:

    // Constructor
	Receiver(double _x, double _y, double _z, double _rad, int _id, double r_mol)
	{
		x = _x;
		y = _y;
		z = _z;

		radius = _rad;

        // Collision radius with a molecule.
		collradius = radius + r_mol;

		id = _id;

		numFalseMolsReceived = 0;
		numMolsReceived = 0;
	}

    // Getter and Setter Functions
	double getX()
	{
		return x;
	}

	double getY()
	{
		return y;
	}

	double getZ()
	{
		return z;
	}

	double getRadius()
	{
		return radius;
	}

    double getCollradius()
    {
        return collradius;
    }

	int getID()
	{
		return id;
	}

	int getNumberOfMoleculesReceived()
	{
		return numMolsReceived;
	}

	int getNumberOfFalseMoleculesReceived()
	{
		return numFalseMolsReceived;
	}

    // Increase Number of Molecules Received that were Destined for another Receiver
    void addFalseReceived()
    {
        this->numFalseMolsReceived++;
    }

    // Increase Number of Molecules Received
    void addReceived()
    {
        this->numMolsReceived++;
    }

private:

    // Position of the Receiver in 3D Space
	double x;
	double y;
	double z;

    // Radius and Collision Radius
	double radius;
	double collradius;

    // ID of the Receiver
	int id;

    // Counters to Hold Number of Received Molecules (Correctly and Wrongly)
	int numMolsReceived;
	int numFalseMolsReceived;
};

#endif
