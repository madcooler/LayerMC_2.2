#include "stdafx.h"
#include "Ray.h"


Ray::Ray(void)
{
	direction=Vector(0,0,0);
	stokes=StokesVector(0,0,0,0);
	Energy=0;
	layerNum=0;
}


Ray::~Ray(void)
{
}


Ray::Ray(Vector d, StokesVector sv,int layer)
{
	direction=d;
	stokes=sv;
	Energy=sv.s[0];
	layerNum=layer;
}
