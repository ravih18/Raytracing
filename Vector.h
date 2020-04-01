#pragma once
#include <math.h>
#include <vector>

#include <random>
static std::default_random_engine engine;
static std::uniform_real_distribution<double> uniform(0, 1);

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	const double& operator[](int i) const { return coord[i];  }
	double& operator[](int i) { return coord[i]; }

	// Give the square norm of the vector
	double getN2() {
		return coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2];
	}

	void normalize() {
		double norm = sqrt(getN2());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
	}
	Vector getNormalized() {
		Vector result(*this);
		result.normalize();
		return result;
		}
	Vector& operator+=(const Vector& b) {
		coord[0] += b[0];
		coord[1] += b[1];
		coord[2] += b[2];
		return *this;
	}

private:
	double coord[3];
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(double a, const Vector& b);
Vector operator*(const Vector& b, double a);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector& a, double b);
double dot(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

Vector random_cos(const Vector& Normal); 

class Ray {
public:
	Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
	Vector origin, direction;
};

class Sphere {
public:
	Sphere(const Vector& origin, double rayon, const Vector &color, bool mirror = false, bool transp = false) : O(origin), R(rayon), albedo(color), isMirror(mirror), isTransp(transp) {};
	
	bool intersect(const Ray& d, Vector& Pos, Vector& Normal, double &t) const {
		// solve at²+bt+c=0
		double a = 1;
		double b = 2 * dot(d.direction, d.origin - O);
		double c = (d.origin - O).getN2() - R * R;
		double delta = b * b - 4 * a * c;

		if (delta < 0) return false;
		double t1 = (-b - sqrt(delta)) / (2 * a);
		double t2 = (-b + sqrt(delta)) / (2 * a);
		// t2 is higher than t1

		if (t2 < 0) return false;
		// from here we know there is an intersection
		if (t1 > 0)
			t = t1;
		else
			t = t2;

		Pos = d.origin + t * d.direction;
		Normal = (Pos - O).getNormalized();

		return true;
	}

	Vector O;
	double R;
	Vector albedo;
	bool isMirror;
	bool isTransp;
};

class Scene {
public:
	Scene() {};
	void addSphere(const Sphere& s) { spheres.push_back(s); }

	bool intersect(const Ray& d, Vector& Pos, Vector& Normal, int &sphere_id, double &t_min) const {

		bool has_inter = false;
		t_min = 1E99;

		for (int i = 0; i < spheres.size(); i++) {
			Vector localPos, localNormal;
			double t;
			bool localHasInter = spheres[i].intersect(d,  localPos, localNormal, t);
			if (localHasInter) {
				has_inter = true;
				if (t < t_min) {
					t_min = t;
					Pos = localPos;
					Normal = localNormal;
					sphere_id = i;
				}
			}
		}
		return has_inter;
	}

	std::vector<Sphere> spheres;
	Sphere *light;
	double light_intensity;
};