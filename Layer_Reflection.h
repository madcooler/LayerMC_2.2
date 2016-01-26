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
#include "stdafx.h"
#include <vector>
#include "Polarization.h"
#include "Ray.h"
#include "Fresnel.h"
//#include "windows.h"
#include "omp.h"
#include <iostream>
#include <fstream>
#include "AnisotropicSampler.h"
using namespace std;




class Layer_Reflection
{
public:
	void AddLayer(float n, float k,float rx,float ry,float lx,float ly);
	void AddAnisotropicLayer(float n, float k,float ex,float ey);
	void SetIncidence(float theta, float phi, StokesVector stokes);
	void CollectData(Ray ray);
	void Compute(void);
	void Calculate(void);

	int NumOfLayer,NumOfIteration;
	Vector Z;
	vector<float> N,K,RLx,RLy,RMSx,RMSy;            //IOR of each layer ,  air as first layer
	StokesVector *reflection_data[91][360];
	
	//float *Normal_X[800*800];//,*Normal_Z[800][800],*Normal_Z[800][800];
	
	Ray OriginIncidence;
	void Trace_Incidence(Ray Incidence,vector<Ray> &IncidenceList);
	Layer_Reflection(void);
	~Layer_Reflection(void);
	Vector SampleLocalNormal(int layer);
	Vector SampleAnisotropicNormal(int layer);
	void Refract(Ray rayin, float n1, float n2, Vector norm,bool &TIR,int layerNum,Ray & rayout);
	
	void Reflection(const Ray in, float n1,float k1, float n2,float k2, Vector norm,int layerNum,Ray & rayout);
	void SaveToFile();
	void ReadSettingFile();
	void ReadAnisotropicSettingFile();

	vector<float> Ex,Ey; 
	
};



