#pragma once
#include <math.h>
#include <float.h>
#define M_PI 3.1415926
#define isnan _isnan
class Vector {
public:
	// Vector Public Methods
	Vector() { x = y = z = 0.f; }
	Vector(float xx, float yy, float zz)
		: x(xx), y(yy), z(zz) {
			
	}
	void Print()
	{
		printf("\n x %f y %f z %f \n", x,y,z);
	}

	// The default versions of these are fine for release builds; for debug
	// we define them so that we can add the Assert checks.
	Vector(const Vector &v) {
		
		x = v.x; y = v.y; z = v.z;
	}

	Vector &operator=(const Vector &v) {
		
		x = v.x; y = v.y; z = v.z;
		return *this;
	}

	Vector operator+(const Vector &v) const {
		
		return Vector(x + v.x, y + v.y, z + v.z);
	}

	Vector& operator+=(const Vector &v) {
		
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	Vector operator-(const Vector &v) const {
		
		return Vector(x - v.x, y - v.y, z - v.z);
	}

	Vector& operator-=(const Vector &v) {
		
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}
	Vector operator*(float f) const { return Vector(f*x, f*y, f*z); }

	Vector &operator*=(float f) {
		
		x *= f; y *= f; z *= f;
		return *this;
	}
	Vector operator/(float f) const {
		if(f==0)
		{
			return Vector(0,0,0);
		}

		
		float inv = 1.f / f;
		return Vector(x * inv, y * inv, z * inv);
	}
	bool HasNaNs() const {
		return isnan(x) || isnan(y) || isnan(z);
	}
	Vector &operator/=(float f) {
		
		float inv = 1.f / f;
		x *= inv; y *= inv; z *= inv;
		return *this;
	}
	Vector operator-() const { return Vector(-x, -y, -z); }
	float operator[](int i) const {
		
		return (&x)[i];
	}

	float &operator[](int i) {
		
		return (&x)[i];
	}
	float LengthSquared() const { return x*x + y*y + z*z; }
	float Length() const { return sqrtf(LengthSquared()); }
	

	bool operator==(const Vector &v) const {
		return x == v.x && y == v.y && z == v.z;
	}
	bool operator!=(const Vector &v) const {
		return x != v.x || y != v.y || z != v.z;
	}
	// Vector Public Data
	float x, y, z;
};


// Geometry Inline Functions


inline Vector operator*(float f, const Vector &v) { return v*f; }
inline float Dot(const Vector &v1, const Vector &v2) {
	
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}




inline float AbsDot(const Vector &v1, const Vector &v2) {
	
	return fabsf(Dot(v1, v2));
}


inline Vector Cross(const Vector &v1, const Vector &v2) {
	
	double v1x = v1.x, v1y = v1.y, v1z = v1.z;
	double v2x = v2.x, v2y = v2.y, v2z = v2.z;
	return Vector((v1y * v2z) - (v1z * v2y),
		(v1z * v2x) - (v1x * v2z),
		(v1x * v2y) - (v1y * v2x));
}


inline Vector Normalize(const Vector &v) { return v / v.Length(); }

inline void CoordinateSystem(const Vector &v1, Vector *v2, Vector *v3) {
	if (fabsf(v1.x) > fabsf(v1.y)) {
		float invLen = 1.f / sqrtf(v1.x*v1.x + v1.z*v1.z);
		*v2 = Vector(-v1.z * invLen, 0.f, v1.x * invLen);
	}
	else {
		float invLen = 1.f / sqrtf(v1.y*v1.y + v1.z*v1.z);
		*v2 = Vector(0.f, v1.z * invLen, -v1.y * invLen);
	}
	*v3 = Cross(v1, *v2);
}

inline Vector SphericalDirection(float sintheta,
	float costheta, float phi) {
		return Vector(sintheta * cosf(phi),
			sintheta * sinf(phi),
			costheta);
}


inline Vector SphericalDirection(float sintheta, float costheta,
	float phi, const Vector &x,
	const Vector &y, const Vector &z) {
		return sintheta * cosf(phi) * x +
			sintheta * sinf(phi) * y + costheta * z;
}

inline float SphericalTheta(const Vector &v) {
	Vector x=Normalize(v);
	
	return acosf(abs(x.z));
}


inline float SphericalPhi(const Vector &v) {
	float p = atan2f(v.y, v.x);
	return (p < 0.f) ? p + 2.f*M_PI : p;
}

template <class T>
inline T Sqr(const T& x) {return x*x;}

inline float Clamp(float val, float low, float high) {
	if (val < low) return low;
	else if (val > high) return high;
	else return val;
}

inline float length(Vector v)
{
	return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float CosBetweenVector(const Vector &v1, const Vector &v2)
{
	return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)/(length(v1)*length(v2));
}

