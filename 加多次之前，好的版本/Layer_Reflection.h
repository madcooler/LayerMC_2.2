//--------------------------------------------------------------
//        layer 0 air
//
//        ````````````````````````````````````````
//        layer 1
//
//        ****************************************
//        layer 2
//        
//        #######################################
//        layer 3  substrate ( Opaque material )
//        
// -------------------------------------------------------------
//        So, the variable NumOfLayer in the code is 3
//
//        Source code developed by Wang Chi from HFUT, China
//--------------------------------------------------------------

#pragma once
#include <vector>
#include "Polarization.h"
#include "Ray.h"
#include "windows.h" 
#include "omp.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#define deg M_PI/180;
#define Rad 180/M_PI;


class Layer_Reflection
{
public:
	void AddLayer(float n, float k,float rx,float ry,float lx,float ly);
	void SetIncidence(float theta, float phi, StokesVector stokes);
	void CollectData(Ray ray);
	void Compute(void);
	void Layer_Reflection::Calculate(void);

	int NumOfLayer,NumOfIteration;
	Vector Z;
	vector<float> N,K,RLx,RLy,RMSx,RMSy;            //IOR of each layer ,  air as first layer
	StokesVector *reflection_data[91][360];
	
	//float *Normal_X[800*800];//,*Normal_Z[800][800],*Normal_Z[800][800];
	
	Ray OriginIncidence;
	void Layer_Reflection::Trace_Incidence(Ray Incidence);
	Layer_Reflection(void);
	~Layer_Reflection(void);
	Vector SampleLocalNormal(int layer);
	Ray Refract(Ray rayin, float n1, float n2, Vector norm,bool &TIR,int layerNum);
	
	Ray Layer_Reflection::Reflection(const Ray in, float n1,float k1, float n2,float k2, Vector norm,int layerNum);
	void SaveToFile();
	void ReadSettingFile();
	
};


inline void Fresnel_Refraction(float theta,float n1,float k1, float n2,float k2,float &phase,float &Ts,float &Tp,float &TsTp,float &factor) 
{
	float nab=n2/n1;
	float theta_t=asin(n1*sin(theta)/n2);
	Tp=(2*nab*cos(theta)/(cos(theta)+nab*cos(theta_t)))*(2*nab*cos(theta)/(cos(theta)+nab*cos(theta_t)));
	Ts=(2*nab*cos(theta)/(nab*cos(theta)+cos(theta_t)))*(2*nab*cos(theta)/(nab*cos(theta)+cos(theta_t)));
	phase=0;
	TsTp=(4*nab*nab*cos(theta)*cos(theta))/((cos(theta)+nab*cos(theta_t))*(nab*cos(theta)+cos(theta_t)));

	float nba=n2/n1;
	factor=nba*nba*nba*(cos(theta_t)/cos(theta));
}

inline float Fresnel_Reflection(float theta,float n1,float k1, float n2,float k2,float &rs,float &rp,float &phaseS,float &phaseP,float &Fs,float &Fp) 
{
	float temp=(n2*n2-k2*k2-n1*n1*sin(theta)*sin(theta));
	float p2 = (sqrt(temp*temp+4*n2*n2*k2*k2)+temp)/2;
	float q2 = (sqrt(temp*temp+4*n2*n2*k2*k2)-temp)/2;
	float p=sqrt(p2),q=sqrt(q2);
	
	Fs=((n1*cos(theta)-p)*(n1*cos(theta)-p)+q2)/((n1*cos(theta)+p)*(n1*cos(theta)+p)+q2);
	Fp=Fs*((p-n1*sin(theta)*tan(theta))*(p-n1*sin(theta)*tan(theta)) +q2)
		/((p+n1*sin(theta)*tan(theta))*(p+n1*sin(theta)*tan(theta)) +q2);

	if (n1>n2&k2==0)
	{
		phaseS=180*deg;
		phaseP=360*deg;
	}
	float ThetaBrewster=atanf(n2/n1);
	if (n1<n2&&k2==0)
	{
		phaseS=180;
		if (theta>ThetaBrewster)
			phaseP=180;
		else phaseP=0;
	}
	if (k2!=0)
	{
		phaseS=atanf(2*n1*q*cos(theta)/(n1*n1*cos(theta)*cos(theta)-p2-q2));
		phaseP=atanf((-2*n1*q*cos(theta)*(p2+q2-n1*n1*sin(theta)*sin(theta)))/((n2*n2+k2*k2)*(n2*n2+k2*k2)*cos(theta)*cos(theta)-n1*n1*(p2+q2)));
	}
	

	rs=sqrt(abs(Fs));
	rp=sqrt(abs(Fp));

	return (Fs+Fp)/2;
}

inline void Fresnel_Reflection(float theta,float n1,float k1, float n2,float k2,STD_COMPLEX &rs,STD_COMPLEX &rp) 
{
	STD_COMPLEX N2(n2,k2),N22;
	STD_COMPLEX N1(n1,k1);
	N22=N2*N2;
	STD_COMPLEX term1,term12,term2,costheta;
	costheta=STD_COMPLEX(cos(theta),0);
	term1=sqrt(N22-STD_COMPLEX(sin(theta)*sin(theta),0));
	term12=N22-STD_COMPLEX(sin(theta)*sin(theta),0);
	term2=sqrt(N22-STD_COMPLEX(n1*n1*sin(theta)*sin(theta),0));
	rs=(costheta-term1)/(costheta+term1);
	rp=(term12*costheta-term2*N1)/(term12*costheta+term2*N1);
	//return (RS+RP)/2;
}



inline double AverageRandom(double min, double max)
{
	int minInteger = (int)(min*10000);
	int maxInteger = (int)(max*10000);
	int randInteger = rand()*rand();
	int diffInteger = maxInteger - minInteger;
	int resultInteger = randInteger % diffInteger + minInteger;
	return resultInteger/10000.0;
}

inline double GaussianNoise(double variance,double mean)
{ 
	
	double Temp;
	double u1,u2;
	u1=AverageRandom(0,1);
	u2=AverageRandom(0,1);
	Temp=sqrt(-2*(log(u1)))*sin(2*M_PI*(u2));
	Temp=variance*Temp+mean;	
	return Temp;
}

