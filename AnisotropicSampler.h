#include "Vector.h"
#include <cstdlib>


inline double RandomNum()
{
	return rand()/(double)RAND_MAX;
}

inline double GaussianNoise(double variance,double mean)
{ 
	
	double Temp;
	double u1,u2;
	u1=RandomNum();
	u2=RandomNum();
	Temp=sqrt(-2*(log(u1)))*sin(2*M_PI*(u2));
	Temp=variance*Temp+mean;	
	return Temp;
}
inline void sampleFirstQuadrant(float ex,float ey,float u1, float u2, float *phi, float *costheta) 
{
		if (ex == ey)
			*phi = M_PI * u1 * 0.5f;
		else
			*phi = atanf(sqrtf((ex+1.f) / (ey+1.f)) *
                     tanf(M_PI * u1 * 0.5f));
		float cosphi = cosf(*phi), sinphi = sinf(*phi);
		*costheta = powf(u2, 1.f/(ex * cosphi * cosphi +
                              ey * sinphi * sinphi + 1));
}

inline void SampleAnisotropic(float ex,float ey, Vector &wh)
{
	float u1;
	float u2;
	u1 = RandomNum();
	u2 = RandomNum();

	float phi, costheta;

    if (u1 < .25f) 
	{
		sampleFirstQuadrant(ex,ey,4.f * u1, u2, &phi, &costheta);
    } 
	else if (u1 < .5f) 
	{
        u1 = 4.f * (.5f - u1);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi = M_PI - phi;
    } 
	else if (u1 < .75f) 
	{
        u1 = 4.f * (u1 - .5f);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi += M_PI;
    } 
	else {
        u1 = 4.f * (1.f - u1);
        sampleFirstQuadrant(ex,ey,u1, u2, &phi, &costheta);
        phi = 2.f * M_PI - phi;
    }
    float sintheta = sqrtf(max(0.f, 1.f - costheta*costheta));
    wh = SphericalDirection(sintheta, costheta, phi);
}




