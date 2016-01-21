#include "stdafx.h"
#include "Polarization.h"

Mueller::
	Mueller(const JonesMat& jj)
{
	//
	// from
	// Bohren and Huffman
	// _Absorption_and_Scattering_of_Light_by_Small_Particles_,
	// (Wiley, New York, 1983). 
	//
	//init();
	double norm0 = norm(jj.j[0]);
	double norm1 = norm(jj.j[1]);
	double norm2 = norm(jj.j[2]);
	double norm3 = norm(jj.j[3]);

	m[0][0] = 0.5*(norm0+norm1+norm2+norm3);
	m[0][1] = 0.5*(norm1-norm0+norm3-norm2);
	m[0][2] = real(jj.j[1]*conj(jj.j[2])+jj.j[0]*conj(jj.j[3]));
	m[0][3] = imag(jj.j[1]*conj(jj.j[2])-jj.j[0]*conj(jj.j[3]));
	m[1][0] = 0.5*(norm1-norm0-norm3+norm2);
	m[1][1] = 0.5*(norm1+norm0-norm3-norm2);
	m[1][2] = real(jj.j[1]*conj(jj.j[2])-jj.j[0]*conj(jj.j[3]));
	m[1][3] = imag(jj.j[1]*conj(jj.j[2])+jj.j[0]*conj(jj.j[3]));
	m[2][0] = real(jj.j[1]*conj(jj.j[3])+jj.j[0]*conj(jj.j[2]));
	m[2][1] = real(jj.j[1]*conj(jj.j[3])-jj.j[0]*conj(jj.j[2]));
	m[2][2] = real(jj.j[0]*conj(jj.j[1])+jj.j[2]*conj(jj.j[3]));
	m[2][3] = imag(jj.j[1]*conj(jj.j[0])+jj.j[3]*conj(jj.j[2]));
	m[3][0] = imag(jj.j[3]*conj(jj.j[1])+jj.j[0]*conj(jj.j[2]));
	m[3][1] = imag(jj.j[3]*conj(jj.j[1])-jj.j[0]*conj(jj.j[2]));
	m[3][2] = imag(jj.j[0]*conj(jj.j[1])-jj.j[2]*conj(jj.j[3]));
	m[3][3] = real(jj.j[0]*conj(jj.j[1])-jj.j[2]*conj(jj.j[3]));
}

Mueller::
	Mueller(const StokesVector& s1,const StokesVector& s2,const StokesVector& s3, const StokesVector& s4)
{
	for (int i=0;i<4;++i) {
		m[0][i]=s1[i];
		m[1][i]=s2[i];
		m[2][i]=s3[i];
		m[3][i]=s4[i];
	}
}

//
// Routine to change Mueller matrix to one for a parity conversion...
//
Mueller 
	Mueller::
	parity() const
{
	Mueller mm = *this;
	mm[0][2] = -mm[0][2];
	mm[0][3] = -mm[0][3];
	mm[1][2] = -mm[1][2];
	mm[1][3] = -mm[1][3];
	mm[2][0] = -mm[2][0];
	mm[2][1] = -mm[2][1];
	mm[3][0] = -mm[3][0];
	mm[3][1] = -mm[3][1];
	return mm;
}

Mueller& 
	Mueller::
	operator=(const Mueller& x) 
{
	for (int i=0;i<4;++i) {
		for (int j=0;j<4;++j) {
			m[i][j]=x.m[i][j];
		}
	}
	return *this;
}

//
// Define multiplication between two Mueller matrices
//
Mueller 
	Mueller::
	operator*(const Mueller& matrix) const 
{
	Mueller out;
	int i,j,k;
	for (i=0;i<4;++i) {
		for (j=0;j<4;++j) { 
			out.m[i][j]=0;
			for (k=0;k<4;++k) {
				out.m[i][j]+=m[i][k]*matrix.m[k][j];
			}
		}
	}
	return out;
}

//
// Define multiplication between a Mueller matrix and a constant...
//
Mueller 
	Mueller::
	operator*(double d) const
{
	Mueller out;

	int i,j;
	for (i=0;i<4;++i) 
		for (j=0;j<4;++j) { 
			out.m[i][j]=m[i][j]*d;
		}
		return out;
}

Mueller 
	Mueller::
	operator/(double d) const
{
	Mueller out;

	int i,j;
	for (i=0;i<4;++i) 
		for (j=0;j<4;++j) { 
			out.m[i][j]=m[i][j]/d;
		}
		return out;
}

//
// Multiplication with a StokesVector ...
//
StokesVector 
	Mueller::
	operator*(const StokesVector& in) const
{
	StokesVector out;
	int i,j;
	for (i=0;i<4;++i) {
		out.s[i]=0;
		for (j=0;j<4;++j) 
			out.s[i]+=m[i][j]*in.s[j];
	}
	out.GetPD();
	return out;
}

//
// Define addition of two Mueller matrices...
//
Mueller 
	Mueller::
	operator+(const Mueller& a) const
{
	Mueller c;

	int i,j;
	for (i=0;i<4;++i) {
		for (j=0;j<4;++j) {
			c.m[i][j] = m[i][j] + a.m[i][j];
		}
	}
	return c;
}

//
// Define subtraction of two Mueller matrices...
//
Mueller 
	Mueller::
	operator-(const Mueller& a) const
{
	Mueller c;

	int i,j;
	for (i=0;i<4;++i) {
		for (j=0;j<4;++j) {
			c.m[i][j] = m[i][j] - a.m[i][j];
		}
	}
	return c;
}

Mueller 
	Mueller::
	operator-() const
{
	Mueller c;

	int i,j;
	for (i=0;i<4;++i) {
		for (j=0;j<4;++j) {
			c.m[i][j] = -m[i][j];
		}
	}
	return c;
}

Mueller Mueller::
	rotate(const double angle) const
{
	Mueller temp(*this);
	return (Mueller(Jones_Rotator(angle))*temp)*
		Mueller(Jones_Rotator(-angle));
}


Mueller 
	Mueller::
	transpose() const
{
	Mueller mm;
	int i,j;
	for (i=0;i<4;++i) 
		for (j=0;j<4;++j) 
			mm[i][j]=m[j][i];
	return mm;
}


Mueller 
	Mueller_Zero()
{
	Mueller m;
	for (int i=0;i<4;++i) 
		for (int j=0;j<4;++j) 
			m[i][j]=0.;
	return m;
}

Mueller 
	MuellerUnit(double transmittance)
{
	Mueller m=Mueller_Zero();
	m[0][0]=m[1][1]=m[2][2]=m[3][3]=transmittance;
	return m;
}

Mueller 
	MuellerDepolarizer(double transmittance,double depolarization)
{
	Mueller m=Mueller_Zero();
	m[0][0]=transmittance;
	m[1][1]=m[2][2]=m[3][3]=transmittance*(1-depolarization);
	return m;
}

Mueller 
	MuellerPartialLinearPolarizer(double tmax, double tmin, double angle)
{
	Mueller m;

	double tpt=0.5*(tmax+tmin);
	double tmt=0.5*(tmax-tmin);
	double tt= sqrt(tmax*tmin);
	double cos2a = cos(2*angle);
	double sin2a = sin(2*angle);

	m[0][0]=tpt;
	m[0][1]=tmt*cos2a;
	m[0][2]=tmt*sin2a;
	m[0][3]=0;
	m[1][0]=m[0][1];
	m[1][1]=tpt*Sqr(cos2a)+tt*Sqr(sin2a); // Correction: TAG 23 JUN 2004
	m[1][2]=(tpt-tt)*cos2a*sin2a;
	m[1][3]=0;
	m[2][0]=m[0][2];
	m[2][1]=m[1][2]; // Correction: TAG 5 DEC 2003
	m[2][2]=tt*Sqr(cos2a)+tpt*Sqr(sin2a);
	m[2][3]=0;
	m[3][0]=m[3][1]=m[3][2]=0;
	m[3][3]=tt;

	return m;
}

double 
	StokesVector ::
	psi() const
{
	return ((JonesVec)((*this).pol_part())).psi();
}

double 
	StokesVector ::
	delta() const
{
	return ((JonesVec)((*this).pol_part())).delta();
}

StokesVector 
	StokesVector ::
	operator*(const Mueller& matrix) const
{
	StokesVector out;
	int i,j;
	for (i=0;i<4;++i) {
		out.s[i]=0;
		for (j=0;j<4;++j) 
			out.s[i]+=s[j]*matrix.m[j][i];
	}
	out.GetPD();
	return out;
}

StokesVector 
	StokesVector ::
	rotate(double angle) const
{
	StokesVector temp(*this);
	return Mueller(Jones_Rotator(angle))*temp;
}

StokesVector 
	StokesVector ::
	pol_part() const
{
	StokesVector ss;
	double dop = DOP();
	ss.s[0] = s[0]*dop;
	ss.s[1] = s[1];
	ss.s[2] = s[2];
	ss.s[3] = s[3];
	ss.GetPD();
	return ss;
}

StokesVector 
	StokesVector ::
	unpol_part() const
{
	StokesVector ss;
	double dop = DOP();
	ss.s[0] = s[0]*(1.-dop);
	ss.s[1] = 0;
	ss.s[2] = 0;
	ss.s[3] = 0;
	ss.GetPD();
	return ss;
}

JonesVec::
	JonesVec(const StokesVector& x)
{
	StokesVector xx = x.pol_part();
	// Next line corrected for Version 4 (TAG 5 MAR 2003)...
	double arg = (xx[3]==0.&&xx[2]==0.) ? 0. : atan2(xx[3],xx[2]);
	double mag = sqrt((xx[0]-xx[1])/2.);
	(*this)[0] = sqrt((xx[0]+xx[1])/2.);
	(*this)[1] = STD_COMPLEX(mag*cos(arg),mag*sin(arg));
}

double 
	JonesVec::
	DOCP() const
{
	return ((StokesVector)(*this)).DOCP();
}

double 
	JonesVec::
	e() const
{
	return ((StokesVector)(*this)).e();
}

double 
	JonesVec::
	DOLP() const
{
	return ((StokesVector)(*this)).DOLP();
}

double 
	JonesVec::
	eta() const
{
	return ((StokesVector)(*this)).eta();
}



JonesVec 
	JonesVec::
	rotate(const double angle) const
{
	return Jones_Rotator(angle)*(*this);
}

 JonesMat::
    JonesMat(const Mueller& mm)
    {
        // 
        // Transformation from Mueller to JonesMat from
        // R.A. Chipman, "Polarimetry" in _Handbook_of_Optics_Vol._II_
        // (McGraw-Hill, New York, 1995).
        // 
        // Modified by TAG to use atan2(,)
        //
        j[0] = sqrt((mm.m[0][0]-mm.m[0][1]-mm.m[1][0]+mm.m[1][1])/2.)*
            exp(STD_COMPLEX(0,1)*atan2((mm.m[3][2]-mm.m[2][3]),(mm.m[2][2]+mm.m[3][3])));
        j[1] = sqrt((mm.m[0][0]+mm.m[0][1]+mm.m[1][0]+mm.m[1][1])/2.);
        j[2] = sqrt((mm.m[0][0]-mm.m[0][1]+mm.m[1][0]-mm.m[1][1])/2.)*
            exp(STD_COMPLEX(0,1)*atan2((-mm.m[0][3]-mm.m[1][3]),(mm.m[0][2]+mm.m[1][2])));
        j[3] = sqrt((mm.m[0][0]+mm.m[0][1]-mm.m[1][0]-mm.m[1][1])/2.)*
            exp(STD_COMPLEX(0,1)*atan2((mm.m[3][0]+mm.m[3][1]),(mm.m[2][0]+mm.m[2][1])));
    }
 
    //******************************************************************************
    //**
    //** The following are the stream i/o functions ...
    //**
    //******************************************************************************

    ostream& operator<<(ostream& os, const JonesMat& j) 
    {
        return os << '(' << j[0] << ',' << j[1] << ',' << j[2] << ',' << j[3] << ')';
    }

    istream& operator>>(istream& is, JonesMat& j) 
    {
        STD_COMPLEX j1,j2,j3,j4;
        char c;
        is >> c;
        if (c=='(') {
            is >> j1;
            is >> c;
            if (c==',') {
                is >> j2;
                is >> c;
                if (c==',') {
                    is >> j3;
                    is >> c;
                    if (c==',') {
                        is >> j4;
                        is >> c;
                        if (c==')') {
                            j = JonesMat(j1,j2,j3,j4);
                            return is;
                        }
                    }
                }
            }
        }
        is.setstate(ios::failbit);
        return is;
    }

    //******************************************************************************
    //**
    //** The following are various optimum polarization modifying elements...
    //**
    //******************************************************************************

    JonesMat 
    Jones_Rotator(double angle)
    {
        double c = cos(angle);
        double s = sin(angle);
        return JonesMat(c,c,s,-s);
    }

    JonesMat 
    JonesLinearRetarder(double phase,double angle) 
    {
        JonesMat out;
        STD_COMPLEX phaser=exp(-STD_COMPLEX(0,1)*phase/2.);
        out[0]=phaser;       // advanced by phase/2
        out[1]=1./phaser;     // delayed by phase/2
        out[2]=out[3]=0;
        return out.rotate(angle);
    }

    JonesMat 
    JonesCircularRetarder(double phase) 
    {
        JonesMat out;
        out[0]=cos(phase/2);
        out[1]=out[0];
        out[2]=sin(phase/2);
        out[3]=-out[2];
        return out;
    }

    JonesMat 
    JonesLinearPolarizer(double angle,double diattenuation)
    {
        JonesMat out(1,sqrt(-(diattenuation-1.)/(diattenuation+1.)),0,0);
        return out.rotate(angle);
    }

    JonesMat
    JonesCircularPolarizer(double diattenuation)
    {
        double e = sqrt((diattenuation-1.)/(diattenuation+1.));
        double x = (1.+e)/2.;
        STD_COMPLEX y = STD_COMPLEX(0.,(1.-e)/2.);
        JonesMat out(x,x,y,-y);
        return out;
    }

    //******************************************************************************
    //**
    //** The following function returns a JonesMat rotated...
    //**
    //******************************************************************************

    JonesMat 
    JonesMat::
    rotate(const double angle) const
    {
        JonesMat rot = Jones_Rotator(angle);
        return rot*((*this)*rot.transpose());
    }

    //******************************************************************************
    //**
    //** The following returns the JonesMat given two Jones eigenvectors and  
    //** their respective eigenvalues (transmittances).
    //**
    //******************************************************************************

    JonesMat 
    JonesGeneralized(const JonesVec& a, const JonesVec& b, 
                     const STD_COMPLEX& ma, const STD_COMPLEX& mb)
    {
        STD_COMPLEX det = a[0]*b[1]-a[1]*b[0];
    
        return JonesMat(a[0]*b[1]*mb-a[1]*b[0]*ma,
                           a[0]*b[1]*ma-a[1]*b[0]*mb,
                           a[0]*b[0]*(mb-ma),
                           a[1]*b[1]*(ma-mb))/det;
    }
