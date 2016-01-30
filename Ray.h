#pragma once
#include "Polarization.h"
#include "Vector.h"

class Ray
{
public:
	Vector direction,p_d,s_d; // propagation and p and s component direction
	StokesVector stokes;
	float Energy;
	int layerNum;
	Ray(void);
	~Ray(void);
	Ray(Vector d, StokesVector sv,int layer);
	
	void SetParameters(Vector d, StokesVector sv,int layer)
	{
		direction=d;
		stokes=sv;
		Energy=sv.s[0];
		layerNum=layer;
	}
	Ray operator = (Ray a)
	{
		direction=a.direction;
		stokes=a.stokes;
		Energy=a.Energy;
		layerNum=a.layerNum;
		return *this;
	}
	bool HasEnergy()
	{
		if(Energy>0.1)
			return true;
		else return false;
	}
	void Print()
	{
		printf("Current Layer %d",layerNum);
		printf("\n direction:");
		direction.Print();
		stokes.Print();
		printf("\n");
	}
};

