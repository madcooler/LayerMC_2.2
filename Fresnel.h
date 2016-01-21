//
//  Fresnel.h
//  MC_Layer
//
//  Created by Administrator on 21/01/16.
//  Copyright © 2016 Administrator. All rights reserved.
//

#ifndef Fresnel_h
#define Fresnel_h


#endif /* Fresnel_h */

#include "Vector.h"

#define deg M_PI/180;
#define Rad 180/M_PI;



inline void GetSpecularRefractionDirection(
            const Vector  incidentDirection,
            const Vector  normalDirection,
            const double  n_i,
            const double  n_t,
            bool & TIR,
            Vector & refractionDirection
            )
{
    
    double eta = n_i/n_t;
    // Make sure the local Normal  is in the same side of incidence
    Vector  localNormal =  normalDirection;
    
    //Propagate Downwards
    double cosi = Dot( incidentDirection ,
                        normalDirection);
    if( cosi > 0 )
        localNormal = - (normalDirection) ;
    
    double idotN = abs(cosi);
    
    double k = 1 - eta*eta*(1 - idotN*idotN); //cos_t * cos_t

	if (k < 0)
	{
        // total internal reflection;
        TIR=true;
//		NSLog(@"TIR");
	}
	else
	{
		TIR=false;
        Vector refraction;
        
        Vector temp1;
        temp1 = eta * (incidentDirection);
        Vector temp2;
        temp2 = (eta*idotN-sqrt(k)) * (localNormal);
        
        refraction = temp1 + temp2 ;
        
        refractionDirection = Normalize(refraction);
        
	}
    
    
    //^^^^^^^ test ^^^^^^^
    /*
    double cost = vec3d_vv_dot( refractionDirection ,
                                normalDirection
                              );
    double sint = sqrt(1-cost*cost);
    double sini = sqrt(1-cosi*cosi);
    NSLog(@" theta_i %f   theta_t  %f",asin(sini),asin(sint));
    NSLog(@" NiSini %f   NtSint  %f",n_i*sini,n_t*sint);
    */
    
    //vec3d_printf(incidentDirection);
    //vec3d_printf(refractionDirection);
    
}

inline void GetSpecularReflectionDirection(
            const Vector incidentDirection,
            const Vector normalDirection,
            const double  n_i,
            const double  n_t,
                  Vector & refractionDirection
            )
{
    Vector  localNormal = normalDirection;
    
    //Propagate Downwards
    double cosi = Dot( incidentDirection ,
                        normalDirection);
    if( cosi > 0 )
        localNormal = - (normalDirection) ;
    
    refractionDirection = Normalize(
                    incidentDirection
                    +
                    2 * (AbsDot(incidentDirection,localNormal)*localNormal)
                    );
}



inline void Fresnel_Refraction(float theta,float n1,float k1, float n2,float k2,float &phase,float &Ts,float &Tp,float &TsTp,float &factor) 
{
	float nab=n1/n2;
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

inline void rotationAngle(
    const   Vector   incomingDirection,
    const   Vector   outgoingDirection,
    const   Vector   surfaceNormal,
    const   Vector   localNormal,
            double & globalTolocalAngle,
            double & localToGlobalAngle
            )
 {
 
    Vector Si,Sr,Ei,Er,temp_in,temp_out;
	Si=Cross(incomingDirection,surfaceNormal);
	Si=Normalize(Si);
	Sr=-Cross(outgoingDirection,surfaceNormal);
	Sr=Normalize(Sr);
	Ei=Cross(incomingDirection,localNormal);
	Ei=Normalize(Ei);
	Er=-Cross(outgoingDirection,localNormal);
	Er=Normalize(Er);
	float cos_rotate_in,cos_rotate_out,sin_rotate_in,sin_rotate_out;
	cos_rotate_in=Clamp(CosBetweenVector(Si,Ei),-1,1);
	cos_rotate_out=Clamp(CosBetweenVector(Sr,Er),-1,1);
	sin_rotate_in=sqrtf(1-cos_rotate_in*cos_rotate_in);
	sin_rotate_out=sqrtf(1-cos_rotate_out*cos_rotate_out);

	//Rotation angle for incident light and refraction light
	globalTolocalAngle=acosf(cos_rotate_in);
	localToGlobalAngle=acosf(cos_rotate_out);

	temp_in=Cross(Si,Ei);
	temp_out=Cross(Er,Sr);   

	//Determine rotating a positive or negtive angle for incident and reflection light
	if (Dot(temp_in,incomingDirection)<0)
	{
		globalTolocalAngle=2*M_PI-(globalTolocalAngle);
	}

	if (Dot(temp_out,outgoingDirection)<0)
	{
		localToGlobalAngle=2*M_PI-(localToGlobalAngle);
	}
 
 
 }

inline bool isNormalReflection(
    const   Vector   incomingDirection,
    const   Vector   outgoingDirection
    )
{
    double zofin  = incomingDirection.z;
    double zofout = outgoingDirection.z;
    
    if (zofin*zofout < 0)
        return true;
    else return false;
}

inline bool isNormalRefraction(
    const   Vector   incomingDirection,
    const   Vector   outgoingDirection
    )
{
    double zofin  = incomingDirection.z;
    double zofout = outgoingDirection.z;
    
    if (zofin*zofout > 0)
        return true;
    else return false;
}