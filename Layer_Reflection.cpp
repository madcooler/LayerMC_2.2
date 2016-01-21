#include "stdafx.h"
#include "Layer_Reflection.h"


Layer_Reflection::Layer_Reflection(void)
{
	N.push_back(1);  // Top layer: AIR
	K.push_back(0);
	NumOfLayer=0;
	NumOfIteration=0;
	
	Z=Vector(0,0,1);

	for (int i=0;i<91;i++)
	{
		for(int j=0;j<360;j++)
			reflection_data[i][j]=new StokesVector();
	}
	
	//ReadSettingFile();
	ReadAnisotropicSettingFile();
	srand(1);
}


Layer_Reflection::~Layer_Reflection(void)
{
}

Vector Layer_Reflection::SampleBlinnNormal(int layer)
{
	int index=layer-1;
	float ex = Ex[index];
	float ey = Ey[index];

	Vector normal;

	SampleAnisotropic(ex,ey,normal);
	
	return normal;

}


//Return the local normal direction according to given number of the ith layer
//using Gaussian sampling
Vector Layer_Reflection::SampleLocalNormal(int layer)
{
	//
	int index=layer-1;
	float x,y,rx,ry,lx,ly;
	rx=RMSx[index];
	ry=RMSy[index];
	lx=RLx[index];
	ly=RLy[index];
	
    // rx is the roughness along the X direction
	x=GaussianNoise(rx,0)-GaussianNoise(rx,0);
	Vector xd(lx,0,x);
    
    // ry is the roughness along the Y direction
	y=GaussianNoise(ry,0)-GaussianNoise(ry,0);
	Vector yd(0,ly,y);

	Vector n=Normalize(Cross(xd,yd));
    
    //make sure the local normal point to the positve direction
	if (n.z<0)
	{
		n=-n;
	}
	if(n.HasNaNs())
	{
		n=SampleLocalNormal(layer);
	}

	//return Normalize(Vector(x,y,z));
	return (n);//return normalized normal direction

	//return Normalize(Vector(1,-1,15)) ;
}

// Add layer parameters including the IOR ,roughness, and 'l' into correspoding stacks
void Layer_Reflection::AddAnisotropicLayer(float n, float k,float ex,float ey)
{
	N.push_back(n);
	K.push_back(k);
	Ex.push_back(ex);
	Ey.push_back(ey);
	
	NumOfLayer++;
}

void Layer_Reflection::AddLayer(float n, float k,float rx,float ry,float lx,float ly)
{
	N.push_back(n);
	K.push_back(k);
	RMSx.push_back(rx);
	RMSy.push_back(ry);
	RLx.push_back(lx);
	RLy.push_back(ly);
	NumOfLayer++;
}

// Set the initial incident light
void Layer_Reflection::SetIncidence(float theta, float phi, StokesVector stokes)
{
	theta=theta*deg;
	phi=phi*deg;
	float sintheta=sin(theta),costheta=cos(theta);
	
	OriginIncidence.direction=Vector(sintheta * cosf(phi),sintheta * sinf(phi),-costheta);
	OriginIncidence.stokes=stokes;
	OriginIncidence.Energy=stokes.s[0];
	OriginIncidence.layerNum=0;
}


void Layer_Reflection::CollectData(Ray ray)
{
	Vector d=Normalize(ray.direction);
	int thetaDeg=0,phiDeg=0;
	
	float x=d.x,y=d.y,z=d.z;
	
	if (z>0)
	{
		
		float costheta=CosBetweenVector(d,Z);
		float temp=acosf(costheta);
		temp=temp*Rad;
		//thetaDeg=(int)temp;
		//int temp1=(int)temp;
		thetaDeg=(int)(temp+0.5);  // >=temp?temp1+1:temp1;
		thetaDeg=Clamp(thetaDeg,0,90);
		

		float absphi=abs(atanf(y/x))*Rad;
		//int temp2=(int)absphi;
		phiDeg=(int)(absphi+0.5);//(2*temp2+1)>=2*absphi?temp2+1:temp2;

		if (x>0&&y>=0)
			phiDeg=phiDeg;
		if (x>0&&y<=0)
			phiDeg=360-phiDeg;
		if (x<0&&y<=0)
			phiDeg=180+phiDeg;
		if (x<0&&y>=0)
			phiDeg=180-phiDeg;

		phiDeg=phiDeg%360;
		phiDeg=Clamp(phiDeg,0,359);

		*reflection_data[thetaDeg][phiDeg]=*reflection_data[thetaDeg][phiDeg]+(ray.stokes);
	}
	
}

void Layer_Reflection::Calculate(void)
{
	#pragma omp parallel for
	for (int i=0;i<NumOfIteration;i++)
	{
		Compute();
	}

	SaveToFile();
}

void Layer_Reflection::Compute(void)
{
	vector<Ray> IncidenceList;
	
	IncidenceList.push_back(OriginIncidence);
	for (int i=0;i<IncidenceList.size();i++)
	{
		Ray incidence=IncidenceList[i];
		Trace_Incidence(incidence,IncidenceList);
	}
	IncidenceList.clear();
}

void Layer_Reflection::Trace_Incidence(Ray Incidence,vector<Ray> &IncidenceList)
{
	int LightLayerNum=Incidence.layerNum;
	/*
	vector<Vector> LocalNormal;
	for (int r=1;r<=NumOfLayer;r++)
	{
		Vector ln=SampleLocalNormal(r); 
		LocalNormal.push_back(ln);
	}
	*/
	
	//Compute the reflection of the first layer
	Ray Top_reflection;
    Reflection(Incidence,N[LightLayerNum],K[LightLayerNum],N[LightLayerNum+1],K[LightLayerNum+1],SampleBlinnNormal(LightLayerNum+1),LightLayerNum,Top_reflection);

	vector<Ray> refraction_ray,reflection_ray,final_reflection;

	//Save the first reflection
    if(isNormalReflection(Incidence.direction, Top_reflection.direction))
        reflection_ray.push_back(Top_reflection);

	bool TIR=false;
	if(LightLayerNum<NumOfLayer-1)
	{
		Ray temp_in=Incidence;
	
		int i=0;
		
		//Calculate all the refractions at each interfaces when the light propagates downwards
		//And store them in "refraction_ray"
        Ray temp_out = Incidence;
		for (i=LightLayerNum;i < NumOfLayer-1;i++)
		{
            temp_in = temp_out;
			Refract(temp_in,N[i],N[i+1],SampleBlinnNormal(i+1),TIR,temp_in.layerNum+1,temp_out);
			if (TIR==false && isNormalRefraction(temp_in.direction, temp_out.direction))
				refraction_ray.push_back(temp_out);
			else break;
		}

		//Calculate all the reflections caused by the refractions (Stored in "refraction_ray") at each interfaces 
		//(except the reflections at the first layer, because it has been calculated above)
		//when the light propagates downwards

		for (int x=0;x<refraction_ray.size();x++ )
		{
			Ray inc=refraction_ray[x];
            Ray temp;
			Reflection(inc,N[inc.layerNum],K[inc.layerNum],N[inc.layerNum+1],K[inc.layerNum+1],SampleBlinnNormal(inc.layerNum+1),inc.layerNum,temp);
			if (isNormalReflection(inc.direction, temp.direction))
				reflection_ray.push_back(temp);
		}
	}

	//Trace the reflections (stored in "reflection_ray") upwards until they leave the top surface
	//Set the new reflections of the reflections (stored in "reflection_ray") as new incident light
	//into the corresponding layer
	for (int i=0;i<reflection_ray.size();i++)
	{
		Ray reftemp=reflection_ray[i];  // Fetch the rays stored in the "reflection_ray"
		int reflayer=reftemp.layerNum;
		TIR=false;
        Ray NewRefraction = reftemp;
		for (int j=reflayer;j>0;j--)
		{
			Vector norm=SampleBlinnNormal(j);
            
			//Calculate the reflections of the ray (reftemp)
			Ray NewReflection;
            reftemp = NewRefraction;
            Reflection(reftemp,N[j],K[j],N[j-1],K[j-1],norm,j,NewReflection);
            
            if(NewReflection.HasEnergy() && isNormalReflection(reftemp.direction, NewReflection.direction))
				IncidenceList.push_back(NewReflection); //Put reflections of the ray(reftemp) into incident list (equals setting it as new incident light)

			//Calculate the refractions of the ray
			
            
            Refract(reftemp,N[j],N[j-1],norm,TIR,reftemp.layerNum-1,NewRefraction);

			
			if (TIR==true || !isNormalRefraction(reftemp.direction, NewRefraction.direction))    
			{
				break;   // if there is Total Internal Reflection happens, just break and stop tracing this ray
			}
		}
		if (TIR==false && isNormalRefraction(reftemp.direction, NewRefraction.direction))
		{
			//cout<<reftemp.Energy<<"   "<<reftemp.direction.z<<endl;

			//Store the rays that survive and leave the surface
			final_reflection.push_back(reftemp); 
		}
		
		//CollectData(reftemp);
	}	

	for (int i=0;i<final_reflection.size();i++)
	{
		CollectData(final_reflection[i]); //Collect the final rays in the upper hemisphere.
	}
	
	//LocalNormal.clear();
	//CollectData(Top_reflection);
}

void Layer_Reflection::Refract(Ray rayin, float n1, float n2, Vector norm,bool &TIR,int layerNum,Ray & rayout)
{
	//Compute refraction direction

	float eta=n1/n2;
	Vector in=Normalize(rayin.direction);
	Vector outdirection;
    
    GetSpecularRefractionDirection(in, norm, n1, n2, TIR, outdirection);
    
    float cosi = Dot(in,norm);
    float IdotN = abs(Dot(in,norm));
	float k = 1 - eta*eta*(1 - IdotN*IdotN);
//
//    Vector normFromTheSameSide = norm;
//    if (cosi > 0)
//    {
//        normFromTheSameSide = - norm;
//    }
	
//
//	if (k<0)
//	{
//		TIR=true;           // total internal reflection;
//		return Ray();
//	} 
//	else
//	{
//		TIR=false;
//		outdirection= eta*in + (eta*IdotN - sqrt(k))*normFromTheSameSide;
//	}
//    
//    outdirection=Normalize(outdirection);

	//Compute refraction stokes
	double Rotate_angle_ref,Rotate_angle_in,cost;
	cosi=abs(IdotN);

	//float sint=eta*eta*(sqrtf(1-cosi*cosi));
	cost=sqrtf(k);

	float theta_i,theta_r;  //local incident and reflection angle
	theta_i=acosf(cosi);
	theta_r=acosf(cost);
	
    rotationAngle(in, outdirection, Z, norm, Rotate_angle_in, Rotate_angle_ref);
	//Calculate reference frame rotation

//	Vector Si,Sr,Ei,Er,temp_in,temp_out;
//	Si=Cross(in,Z);
//	Si=Si/length(Si);
//	Sr=-Cross(outdirection,Z);
//	Sr=Sr/length(Sr);
//	Ei=Cross(in,norm);
//	Ei=Ei/length(Ei);
//	Er=-Cross(outdirection,norm);
//	Er=Er/length(Ei);
//	float cos_rotate_in,cos_rotate_out,sin_rotate_in,sin_rotate_out;
//	cos_rotate_in=Clamp(CosBetweenVector(Si,Ei),-1,1);
//	cos_rotate_out=Clamp(CosBetweenVector(Sr,Er),-1,1);
//	sin_rotate_in=sqrtf(1-cos_rotate_in*cos_rotate_in);
//	sin_rotate_out=sqrtf(1-cos_rotate_out*cos_rotate_out);
//
//	//Rotation angle for incident light and refraction light
//	Rotate_angle_in=acosf(cos_rotate_in);  
//	Rotate_angle_ref=acosf(cos_rotate_out);
//
//	temp_in=Cross(Si,Ei);
//	temp_out=Cross(Er,Sr);
//
//	//Determine rotating a positive or negtive angle for incident and reflection light
//	if (Dot(temp_in,in)<0)
//	{
//		Rotate_angle_in=2*M_PI-Rotate_angle_in;
//	}
//
//	if (Dot(temp_out,outdirection)>0)
//	{
//		Rotate_angle_ref=2*M_PI-Rotate_angle_ref;
//	}

	//Rotate_angle_in=acosf((Dot(norm,Z)-(Dot(in,Z))*Dot(in,norm))/(sin(theta_i)*sin(localtheta)));
	//Rotate_angle_ref=acosf((Dot(norm,Z)-(Dot(outdirection,Z))*Dot(outdirection,norm))/(sin(theta_r)*sin(localtheta)));

	//Compute ts tp
	

	float Ts,Tp,TsTp,phase,factor;
	Fresnel_Refraction(theta_i,n1,0,n2,0,phase,Ts,Tp,TsTp,factor);
	
	
	float A,B,C,S;
	A=(Ts+Tp)/2;
	B=(Ts-Tp)/2;
	C=TsTp;
	S=0;

	JonesMat R1=Jones_Rotator(Rotate_angle_in),R2=Jones_Rotator(Rotate_angle_ref);

	Mueller mu,mR1(R1),mR2(R2);
	mu.SetValue(0,0,A);
	mu.SetValue(1,1,A);
	mu.SetValue(0,1,B);
	mu.SetValue(1,0,B);
	mu.SetValue(2,2,C);
	mu.SetValue(3,3,C);
	mu.SetValue(2,3,S);
	mu.SetValue(3,2,-S);
	mu=mu*factor;
	mu=mR2*mu;
	mu=mu*mR1;
	StokesVector outstokes=mu*rayin.stokes;
	//float energy=outstokes[0];

	Ray RefractionRay(outdirection,outstokes,layerNum);

	rayout = RefractionRay;
}

void Layer_Reflection::Reflection(const Ray in, float eta_i,float k1, float eta_t,float k2, Vector norm,int layerNum,Ray & rayout)
{
	// 从下往上反射时候，法向量应该和入射光在同侧
	Vector outdirection;
//    =Normalize(in.direction + 2*(AbsDot(in.direction,norm)*norm));   //reflection direction
    GetSpecularReflectionDirection(in.direction, norm, eta_i, eta_t, outdirection);
	
	//k2
	bool IsMetal=(k2!=0);
	
	double cosi,cost,localtheta,Rotate_angle_ref,Rotate_angle_in;
	
	cosi=abs(Dot(in.direction,norm));
	float sint=eta_i/eta_t*(sqrtf(1-cosi*cosi));
	cost=sqrtf(1-sint*sint);

	float theta_i,theta_r;
	theta_i=acosf(cosi);
	theta_r=acosf(cost);

	localtheta=acos(Dot(norm,(in.direction)));  //local incident angle
    rotationAngle(in.direction, outdirection, Z, norm, Rotate_angle_in, Rotate_angle_ref);

//	Vector Si,Sr,Ei,Er,temp_in,temp_out;
//	Si=Cross(in.direction,Z);
//	Si=Normalize(Si);
//	Sr=-Cross(outdirection,Z);
//	Sr=Normalize(Sr);
//	Ei=Cross(in.direction,norm);
//	Ei=Normalize(Ei);
//	Er=-Cross(outdirection,norm);
//	Er=Normalize(Er);
//	float cos_rotate_in,cos_rotate_out,sin_rotate_in,sin_rotate_out;
//	cos_rotate_in=Clamp(CosBetweenVector(Si,Ei),-1,1);
//	cos_rotate_out=Clamp(CosBetweenVector(Sr,Er),-1,1);
//	sin_rotate_in=sqrtf(1-cos_rotate_in*cos_rotate_in);
//	sin_rotate_out=sqrtf(1-cos_rotate_out*cos_rotate_out);
//
//	//Rotation angle for incident light and refraction light
//	Rotate_angle_in=acosf(cos_rotate_in);
//	Rotate_angle_ref=acosf(cos_rotate_out);
//
//	temp_in=Cross(Si,Ei);
//	temp_out=Cross(Er,Sr);   
//
//	//Determine rotating a positive or negtive angle for incident and reflection light
//	if (Dot(temp_in,in.direction)<0)
//	{
//		Rotate_angle_in=2*M_PI-Rotate_angle_in;
//	}
//
//	if (Dot(temp_out,outdirection)>0)
//	{
//		Rotate_angle_ref=2*M_PI-Rotate_angle_ref;
//	}

	float rs,rp,phaseS,phaseP,Fs,Fp;
	Fresnel_Reflection(theta_i,eta_i,k1,eta_t,k2,rs,rp,phaseS,phaseP,Fs,Fp);

	float A,B,C,S;
	A=(Fs+Fp)/2;
	B=(Fs-Fp)/2;
	C=cos(phaseS-phaseP)*sqrt(Fs*Fp);
	S=sin(phaseS-phaseP)*sqrt(Fs*Fp);

	JonesMat R1=Jones_Rotator(Rotate_angle_in),R2=Jones_Rotator(Rotate_angle_ref);

	Mueller mu,mR1(R1),mR2(R2);
	mu.SetValue(0,0,A);
	mu.SetValue(1,1,A);
	mu.SetValue(0,1,B);
	mu.SetValue(1,0,B);
	mu.SetValue(2,2,C);
	mu.SetValue(3,3,C);
	mu.SetValue(2,3,S);
	mu.SetValue(3,2,-S);
	mu=mR2*mu;
	mu=mu*mR1;

	StokesVector outstokes=mu*in.stokes;
	

	Ray ReflectionRay(outdirection,outstokes,layerNum);

	rayout = ReflectionRay;
}

void Layer_Reflection::SaveToFile()
{
	char I[]="I.txt";
	char Q[]="Q.txt";
	char U[]="U.txt";
	char V[]="V.txt";
	char DOP[]="DOP.txt";
	char AOP[]="AOP.txt";

	ofstream Ifs(I),Qfs(Q),Ufs(U),Vfs(V),DOPfs(DOP),AOPfs(AOP);

	for (int i=0;i<91;i++)
	{
		for(int j=0;j<360;j++)
		{
			StokesVector s=*reflection_data[i][j];
			Ifs<<s.I()<<" ";
			Qfs<<s.Q()<<" ";
			Ufs<<s.U()<<" ";
			Vfs<<s.V()<<" ";
			DOPfs<<s.DOP()<<" ";
			AOPfs<<s.eta()<<" ";
		}
		Ifs<<endl;
		Qfs<<endl;
		Ufs<<endl;
		Vfs<<endl;
		DOPfs<<endl;
		AOPfs<<endl;
	}

}

void Layer_Reflection::ReadSettingFile()
{
	float incident,I,Q,U,V;
		
	ifstream x("setting.txt");
	if (!x.is_open())
	{
		cout<<"Couldn't find SETTING files"<<endl;
	}
	
	x>>incident;

	x>>I>>Q>>U>>V;

	this->SetIncidence(incident,0,StokesVector(I,Q,U,V));

	x>>NumOfIteration;
	
	cout<<"Incident Angle : "<<incident<<endl;
	cout<<"Incident Stokes : "<<I<<" "<<Q<<" "<<U<<" "<<V<<endl;
	cout<<"Number Of Iteration : "<<NumOfIteration<<endl;

	int c;
	x>>c;
	
	cout<<"Number of Layer : "<<c<<endl;
	cout<<"Layer 0 (AIR)  IOR: ( 1, 0 )"<<endl;
	for (int j=0;j<c;j++)
	{
		float n,k,ly,lx,rx,ry;
		x>>n>>k>>rx>>ry>>lx>>ly;
		
		
		this->AddLayer(n,k,rx,ry,lx,ly);
		cout<<"Layer "<<j+1<< " IOR: ( "<<n<<", "<<k<<" )"<<endl;
	}
	
	x.close();
}

void Layer_Reflection::ReadAnisotropicSettingFile()
{
	float incident,I,Q,U,V;
		
	ifstream x("settingani.txt");
	if (!x.is_open())
	{
		cout<<"Couldn't find SETTING files"<<endl;
	}
	
	x>>incident;

	x>>I>>Q>>U>>V;

	this->SetIncidence(incident,0,StokesVector(I,Q,U,V));

	x>>NumOfIteration;
	
	cout<<"Incident Angle : "<<incident<<endl;
	cout<<"Incident Stokes : "<<I<<" "<<Q<<" "<<U<<" "<<V<<endl;
	cout<<"Number Of Iteration : "<<NumOfIteration<<endl;

	int c;
	x>>c;
	
	cout<<"Number of Layer : "<<c<<endl;
	cout<<"Layer 0 (AIR)  IOR: ( 1, 0 )"<<endl;
	for (int j=0;j<c;j++)
	{
		float n,k,ex,ey;
		x>>n>>k>>ex>>ey;
		
		
		this->AddAnisotropicLayer(n,k,ex,ey);
		cout<<"Layer "<<j+1<< " IOR: ( "<<n<<", "<<k<<" )"<<endl;
	}
	
	x.close();
}





