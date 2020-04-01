// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <vector>
#include <algorithm>
#include "Vector.h"

#define M_PI 3.1415926535897932


inline double sqr(double x) { return x * x; }

void save_image(const char* filename, const unsigned char* tableau, int w, int h) { // (0,0) is top-left corner

	FILE* f;

	int filesize = 54 + 3 * w * h;

	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);
	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	f = fopen(filename, "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	unsigned char* row = new unsigned char[w * 3];
	for (int i = 0; i < h; i++)
	{
		for (int j = 0; j < w; j++) {
			row[j * 3] = tableau[(w * (h - i - 1) * 3) + j * 3 + 2];
			row[j * 3 + 1] = tableau[(w * (h - i - 1) * 3) + j * 3 + 1];
			row[j * 3 + 2] = tableau[(w * (h - i - 1) * 3) + j * 3];
		}
		fwrite(row, 3, w, f);
		fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
	}
	fclose(f);
	delete[] row;
}

Vector getPixel(const Ray& r, const Scene& s, int n_max) {
	// n_max : maximum number of rebounds 

	if (n_max == 0) return Vector(0., 0., 0.);

	Vector Pos, Normal;
	int sphere_id;
	double t;
	bool has_inter = s.intersect(r, Pos, Normal, sphere_id, t);

	Vector pixel(0., 0., 0.);
	if (has_inter) {

		/*if (sphere_id == 0) {
			return s.light->albedo * s.light_intensity / (4 * M_PI * sqr(s.light->R));
		};*/

		if (s.spheres[sphere_id].isMirror) {
			Vector d_mirror = r.direction - 2 * dot(Normal, r.direction) * Normal;
			Ray r_mirror(Pos + 0.001 * Normal, d_mirror);

			pixel = getPixel(r_mirror, s, n_max - 1);
		}
		else if (s.spheres[sphere_id].isTransp) {
			double n1 = 1;
			double n2 = 1.3;
			Vector transpNormal(Normal);
			if (dot(r.direction, Normal) > 0) { //we have to test if we go in the sphere or out
				n1 = 1.3;
				n2 = 1;
				transpNormal = -1.0 * Normal;
			}
			double rad = 1 - sqr(n1 / n2) * (1 - sqr(dot(transpNormal, r.direction)));
			if (rad > 0) {
				Vector d_refract = (n1 / n2) * (r.direction - dot(r.direction, transpNormal) * transpNormal) - transpNormal * sqrt(rad);
				Ray r_refract(Pos + 0.001 * transpNormal, d_refract);

				pixel = getPixel(r_refract, s, n_max - 1);
			}
		}
		else {
			// direct lighting
			//Ray ray_light(Pos + 0.01 * Normal, (s.light_pos - Pos).getNormalized());
			//Vector P_light, N_light;
			//int light_id;

			//double t_light; // distance between intersection and light
			//bool has_inter_light = s.intersect(ray_light, P_light, N_light, light_id, t_light);
			//double d_light2 = (s.light_pos - Pos).getN2(); // distance between point and light
			//if (!(has_inter_light && t_light * t_light < d_light2 * 0.99)) {
			//	double pixel_intensity = s.light_intensity * dot((s.light_pos - Pos).getNormalized(), Normal) / d_light2;
			//	pixel = s.spheres[sphere_id].albedo / M_PI * pixel_intensity;
			//}
			Vector axeOP = (Pos - s.light->O).getNormalized();
			Vector r_dir = random_cos(axeOP);
			Vector random_point = r_dir * s.light->R + s.light->O;
			Vector wi = (random_point - Pos).getNormalized();
			double d_light2 = (random_point - Pos).getN2();
			Vector Normal2 = r_dir;

			Ray ray_light(Pos + 0.01 * Normal, wi);
			Vector P_light, N_light;
			int light_id;

			double t_light; // distance between intersection and light
			bool has_inter_light = s.intersect(ray_light, P_light, N_light, light_id, t_light);

			if (has_inter_light && t_light * t_light < d_light2 * 0.99) {
				pixel = Vector(0., 0., 0.);
			}
			else {
				pixel = s.light_intensity / (4 * M_PI * d_light2) * std::max(0., dot(Normal, wi)) * dot(Normal2, wi) / dot(axeOP, r_dir) * s.spheres[sphere_id].albedo;
			}

			//indirect lighting

			Vector rand_dir = random_cos(Normal);
			Ray r_rand(Pos + 0.01 * Normal, rand_dir);

			pixel += getPixel(r_rand, s, n_max - 1) * s.spheres[sphere_id].albedo;
		}
	}
	return pixel;
}

		

int main()
{
	int W = 1024;
	int H = 1024;
	double fov = 60 * M_PI / 180;
	const int nrays = 8;

	Sphere lum(Vector(20, 70, -30), 5, Vector(1.,1.,1.));
	Sphere sphere1(Vector(0, 0, -70), 20, Vector(1, 0, 1));
	Sphere grnd(Vector(0, -2000 - 30, 0), 2000, Vector(1, 1, 1));
	Sphere roof(Vector(0, 2000 + 80, 0), 2000, Vector(1, 0, 0));
	Sphere right(Vector(2000 + 60, 0, 0), 2000, Vector(0, 1, 0));
	Sphere left(Vector(-2000-60, 0, 0), 2000, Vector(0, 0, 1));
	Sphere back(Vector(0, 0, -2000-150), 2000, Vector(0, 1, 1));

	Scene s;
	s.addSphere(lum);

	s.addSphere(sphere1);
	s.addSphere(grnd);
	s.addSphere(roof);
	s.addSphere(right);
	s.addSphere(left);
	s.addSphere(back);

	s.light = &lum;
	s.light_intensity = 1000000000;

	std::vector<unsigned char> img(W * H * 3);

#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			Vector pixel(0., 0., 0.);
			for (int k = 0; k < nrays; k++) {
				// Box-Muller Method for anti aliasing
				double r1 = uniform(engine);
				double r2 = uniform(engine);
				double R = sqrt(-2 * log(r1));
				double dx = R * cos(2 * M_PI * r2);
				double dy = R * sin(2 * M_PI * r2);

				Vector direction(j - W / 2 + 0.5 + dx, i - H / 2 + 0.5 + dy, -W / (2 * tan(fov / 2)));
				direction.normalize();

				Ray r(Vector(0, 0, 0), direction);

				pixel += getPixel(r, s, 5) / nrays;
			}

			img[((H-i-1) * W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(pixel[0], 1 / 2.2))); //red
			img[((H-i-1) * W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(pixel[1], 1 / 2.2)));; //green 
			img[((H-i-1) * W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(pixel[2], 1 / 2.2)));; //blue
		}
	}


	save_image("out.bmp", &img[0], W, H);

	return 0;
}



