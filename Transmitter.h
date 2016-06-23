#ifndef INC_TRANSMITTER_H
#define INC_TRANSMITTER_H

class Transmitter
{
public:

    // Constructors
	Transmitter(double _x, double _y, double _z, int _id, int _dest, double r_mol)
	{
		x = _x;
		y = _y;
		z = _z; 

		id = _id;
		destination = _dest;

		radius = 10;
		collradius = radius + r_mol;
	}

	Transmitter(double _x, double _y, double _z, double _rad, int _id, int _dest, double r_mol)
	{
		x = _x;
		y = _y;
		z = _z;

		id = _id;
		destination = _dest;

		radius = _rad;
		collradius = radius + r_mol;
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

	int getDestination()
	{
		return destination;
	}

private:

    // Position of the Transmitter in 3D Space
	double x;
	double y;
	double z;

    // Radius and Collision Radius
	double radius;
	double collradius;

    // ID of the Transmitter
	int id;

    // ID of the Receiver that the Transmitter Sends Molecules to.
	int destination;
};

#endif
