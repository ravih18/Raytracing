#pragma once
#include <math.h>
#include "Vector.h"

#define M_PI 3.1415926535897932


// On définit les opérateurs pour les opérations avec les éléments de la classe Vector
Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(double a, const Vector &b) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &b, double a) {
	return Vector(a * b[0], a * b[1], a * b[2]);
}
Vector operator*(const Vector &a, const Vector &b) {
	return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}
Vector operator/(const Vector& a, double b) {
	return Vector(a[0] / b, a[1] / b, a[2] / b);
}

double dot(const Vector& a, const Vector& b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

Vector random_cos(const Vector& Normal) {
	double r1 = uniform(engine);
	double r2 = uniform(engine);
	Vector rand_dir_local(cos(2 * M_PI * r1) * sqrt(1 - r2), sin(2 * M_PI * r1) * sqrt(1 - r2), sqrt(r2));
	Vector random(uniform(engine) - 0.5, uniform(engine) - 0.5, uniform(engine) - 0.5);
	Vector tan1 = cross(Normal, random); tan1.normalize();
	Vector tan2 = cross(tan1, random);

	return   rand_dir_local[2] * Normal + rand_dir_local[0] * tan1 + rand_dir_local[1] * tan2;
}