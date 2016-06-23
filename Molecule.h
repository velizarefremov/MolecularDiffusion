#ifndef INC_MOLECULE_H
#define INC_MOLECULE_H

class Molecule
{
public:

    // Constructors
	Molecule()
	{

	}

	Molecule(double _x, double _y, double _z, int _destid)
	{
		x = _x;
		y = _y;
		z = _z;

		destID = _destid;
	}

    // Getter and Setter Functions
	int getDestination()
	{
		return destID;
	}

	void setDestination(int destid)
	{
		destID = destid;
	}

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

    void setX(double _x)
    {
        this->x = _x;
    }

    void setY(double _y)
    {
        this->y = _y;
    }

    void setZ(double _z)
    {
        this->z = _z;
    }

    // Set Position of the Molecule in 3D Space
    void setPosition(double _x, double _y, double _z)
    {
        this->x = _x;
        this->y = _y;
        this->z = _z;
    }

    // Check whether the molecule is outside the environment boundary.
    // Currently a spherical environment is assumed.
	bool checkBoundary(double r_env)
	{
		if ((x*x + y*y + z*z) > (r_env * r_env))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

private:

    // Position of the Molecule in 3D Space
	double x;
	double y;
	double z;

    // ID of the Destination Receiver.
	int destID;
};

#endif
