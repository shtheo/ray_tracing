#define _CRT_SECURE_NO_WARNINGS 1

#include <vector>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <random>
#include "CAgp9r15.h"
#include <string>
#include <stdio.h>
#include <sstream>

using namespace std;

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define M_PI 3.141592653589793238

std::default_random_engine e;
std::uniform_real_distribution<double> U(0, 1);

class Vector {
public:
	Vector(double x = 0, double y = 0, double z = 0) :x(x), y(y), z(z) {};
	double norm2() { return x*x + y*y + z*z; }
	void normalize() {
		double n = sqrt(norm2());
		x /= n;
		y /= n;
		z /= n;
	}
	double x, y, z;
};

std::ostream &operator<<(std::ostream &os, Vector v) {
	std::stringstream ss;
	ss << "Vector(" << v.x << ", " << v.y << ", " << v.z << ")";
	return os << ss.str();
}

Vector operator+(const Vector& a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}
Vector operator-(const Vector& a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}
Vector operator*(const double& a, const Vector &b) {
	return Vector(a * b.x, a * b.y, a * b.z);
}
Vector operator*(const Vector& a, const double &b) {
	return b * a;
}
Vector operator*(const Vector& a, const Vector& b) {
		return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}
Vector operator/(const Vector& a, const double &b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}
double dot(const Vector& a, const Vector& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}

class Ray {
public:
	Ray(const Vector& C, const Vector& u) : C(C), u(u) {};
	Vector u, C;
};

class Object {
public:
	Object() {};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) = 0;
	Vector albedo;
	bool mirror;
	bool transparent;
	bool light;
};

class Sphere : public Object {
public:
	Sphere(const Vector& O, double R, const Vector& albedo, bool mirror, bool transparent, bool light = false) : 
		O(O), R(R) 
		{
			this->albedo = albedo;
			this->mirror = mirror;
			this->transparent = transparent;
			this->light = light;
		};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		// résolution de t²|u|² + 2t*(u scalaire C-O) + |C-O|² = R²
		color = this->albedo;
		double a = 1;
		double b = 2 * dot(r.u, r.C - O);
		double c = (r.C - O).norm2() - R * R;
		double delta = b * b - 4 * a*c;
		if (delta >= 0) {
			double t1 = (-b - sqrt(delta)) / (2 * a);
			double t2 = (-b + sqrt(delta)) / (2 * a);
			if (t1 > 0) {
				P = t1 * r.u + r.C;
				N = P - O;
				N.normalize();
			}
			else {
				if (t2 > 0) {
					P = t2 * r.u + r.C;
					N = P - O;
					N.normalize();
				}
				else {
					return false;
				}
			}
		}
		return delta >= 0;
	};
	Vector O;
	double R;
};

class Triangle : public Object {
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& albedo, bool mirror, bool transparent, bool light = false) :
		A(A), B(B), C(C)
	{
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
		this->light = light;
	};
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		color = this->albedo;
		N = cross(B - A, C - A);
		N.normalize();
		if (dot(r.u, N) > 0)
		{
			N = -1. * N;
		}
		double t = dot(A - r.C, N) / dot(r.u, N);
		P = t * r.u + r.C;
		double det_matrix = (B - A).norm2()*(C - A).norm2() - pow(dot(B - A, C - A), 2.);
		double beta = (dot(P - A, B - A)*(A - C).norm2() - dot(C - A, B - A)*dot(P - A, C - A)) / det_matrix;
		double gamma = ((B - A).norm2()*dot(P - A, C - A) - dot(B - A, C - A)*dot(P - A, B - A)) / det_matrix;
		double alpha = 1 - beta - gamma;
		if (beta >= 0 && beta <= 1 && gamma >= 0 && gamma <= 1 && alpha >= 0 && alpha <= 1 && t >= 0) {
			return true;
		}
		else {
			return false;
		}
	};
	void get_barycentric(double& alpha, double& beta, double& gamma, const Vector& P) {
		double det_matrix = (B - A).norm2()*(C - A).norm2() - pow(dot(B - A, C - A), 2.);
		beta = (dot(P - A, B - A)*(A - C).norm2() - dot(C - A, B - A)*dot(P - A, C - A)) / det_matrix;
		gamma = ((B - A).norm2()*dot(P - A, C - A) - dot(B - A, C - A)*dot(P - A, B - A)) / det_matrix;
		alpha = 1 - beta - gamma;
	}
	Vector A;
	Vector B;
	Vector C;
};

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class BBox{
public:
	BBox() {};
	BBox(const Vector& bmin, const Vector& bmax) : bmin(bmin), bmax(bmax) {};
	virtual bool intersect(const Ray& r) {
		double tx_1 = (bmin - r.C).x / r.u.x;
		double tx_2 = (bmax - r.C).x / r.u.x;
		double ty_1 = (bmin - r.C).y / r.u.y;
		double ty_2 = (bmax - r.C).y / r.u.y;
		double tz_1 = (bmin - r.C).z / r.u.z;
		double tz_2 = (bmax - r.C).z / r.u.z;
		double tx_min = std::min(tx_1, tx_2);
		double tx_max = std::max(tx_1, tx_2);
		double ty_min = std::min(ty_1, ty_2);
		double ty_max = std::max(ty_1, ty_2);
		double tz_min = std::min(tz_1, tz_2);
		double tz_max = std::max(tz_1, tz_2);
		double min_global = std::max(tx_min, std::max(ty_min, tz_min));
		double max_global = std::min(tx_max, std::min(ty_max, tz_max));
		if (min_global <= max_global) {
			return true;
		}
		else {
			return false;
		}
	}
	double get_plan_separateur(Vector& axe) {
		if ((bmax.x - bmin.x) > (bmax.y - bmin.y) && (bmax.x - bmin.x) > (bmax.z - bmin.z))
		{
			axe = Vector(1, 0, 0);
			return (bmax.x + bmin.x) / 2;
		}
		else {
			if ((bmax.y - bmin.y) > (bmax.x - bmin.x) && (bmax.y - bmin.y) > (bmax.z - bmin.z))
			{
				axe = Vector(0, 1, 0);
				return (bmax.y + bmin.y) / 2;
			}
			else {
				axe = Vector(0, 0, 1);
				return (bmax.z + bmin.z) / 2;
			}
		}
	};
	Vector bmin, bmax;
};

class BVHNode {
public:
	BVHNode() {};
	BVHNode(int debut, int fin, BBox b) : debut(debut), fin(fin), b(b){};
	bool intersect(const Ray& r, std::vector<BVHNode*>& leaves) {
		if (b.intersect(r)) {
			if (!fg || !fd) {

				leaves.push_back(this);
				return true;
			}
			bool found_fg = fg->intersect(r, leaves);
			bool found_fd = fd->intersect(r, leaves);
			return found_fd || found_fg;
		}
		return false;
	}
	BVHNode *fg, *fd;
	BBox b;
	int debut, fin;
};

class Geometry : public Object {
public:
	Geometry(const Vector& offset, double sampling, const Vector& albedo, bool mirror, bool transparent, bool light = false) {
		this->offset = offset;
		this->sampling = sampling;
		this->albedo = albedo;
		this->mirror = mirror;
		this->transparent = transparent;
		this->light = light;
	};

	void readOBJ(const char* obj, bool load_textures) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec.x, &vec.y, &vec.z, &col.x, &col.y, &col.z) == 6) {
					col.x = std::min(1., std::max(0., col.x));
					col.y = std::min(1., std::max(0., col.y));
					col.z = std::min(1., std::max(0., col.z));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				}
				else {
					sscanf(line, "v %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec.x, &vec.y, &vec.z);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec.x, &vec.y);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				}
				else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					}
					else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						}
						else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					}
					else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						}
						else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							}
							else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								}
								else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}

			}

		}
		fclose(f);
		for (int i = 0; i < vertices.size(); i++)
		{
			vertices[i].x = vertices[i].x * sampling + offset.x;
			double new_y = vertices[i].z * sampling + offset.y;
			vertices[i].z = vertices[i].y * sampling + offset.z;
			vertices[i].y = new_y;
		}
		bb = create_box(0, indices.size());
		BVHNode first_node(0, indices.size(), bb);
		geometry_node = first_node;
		cout << "Building the nodes..." << endl;
		build_node(geometry_node);
	}
	void build_node(BVHNode& node) {
		if (node.fin - node.debut > 5) {
			Vector axe;
			double value = node.b.get_plan_separateur(axe);
			int pivot = node.debut;
			for (int i = node.debut; i < node.fin; i++)
			{
				Vector center = (vertices[indices[i].vtxi] + vertices[indices[i].vtxj] + vertices[indices[i].vtxk]) / 3;
				if (dot(center, axe) < value)
				{
					std::swap(indices[i], indices[pivot]);
					pivot++;
				}
			}
			if (pivot != node.fin && pivot != node.debut) {
				BBox fg_box = create_box(node.debut, pivot);
				node.fg = new BVHNode(node.debut, pivot, fg_box);
				build_node(*node.fg);
				BBox fd_box = create_box(pivot, node.fin);
				node.fd = new BVHNode(pivot, node.fin, fd_box);
				build_node(*node.fd);
			}
		}
	}
	BBox create_box(int debut, int fin) {
		BBox b;
		double xmax = vertices[indices[debut].vtxi].x;
		double xmin = vertices[indices[debut].vtxi].x;
		double ymax = vertices[indices[debut].vtxi].y;
		double ymin = vertices[indices[debut].vtxi].y;
		double zmax = vertices[indices[debut].vtxi].z;
		double zmin = vertices[indices[debut].vtxi].z;
		for (int i = debut; i < fin; i++)
		{
			xmax = std::max(vertices[indices[i].vtxi].x, xmax);
			xmax = std::max(vertices[indices[i].vtxj].x, xmax);
			xmax = std::max(vertices[indices[i].vtxk].x, xmax);
			xmin = std::min(vertices[indices[i].vtxi].x, xmin);
			xmin = std::min(vertices[indices[i].vtxj].x, xmin);
			xmin = std::min(vertices[indices[i].vtxk].x, xmin);
			ymax = std::max(vertices[indices[i].vtxi].y, ymax);
			ymax = std::max(vertices[indices[i].vtxj].y, ymax);
			ymax = std::max(vertices[indices[i].vtxk].y, ymax);
			ymin = std::min(vertices[indices[i].vtxi].y, ymin);
			ymin = std::min(vertices[indices[i].vtxj].y, ymin);
			ymin = std::min(vertices[indices[i].vtxk].y, ymin);
			zmax = std::max(vertices[indices[i].vtxi].z, zmax);
			zmax = std::max(vertices[indices[i].vtxj].z, zmax);
			zmax = std::max(vertices[indices[i].vtxk].z, zmax);
			zmin = std::min(vertices[indices[i].vtxi].z, zmin);
			zmin = std::min(vertices[indices[i].vtxj].z, zmin);
			zmin = std::min(vertices[indices[i].vtxk].z, zmin);
		}
		b.bmax = Vector(xmax, ymax, zmax);
		b.bmin = Vector(xmin, ymin, zmin);
		return b;
	}
	virtual bool intersect(const Ray& r, Vector& P, Vector& N, Vector& color) {
		int debut, fin;
		std::vector<BVHNode*> leaves;
		if (!geometry_node.intersect(r, leaves)) {
			return false;
		}
		else {
			bool found = false;
			double dist_min = 1E99;
			int debut, fin;
			for (int j = 0; j < leaves.size(); j++) {
				debut = leaves[j]->debut;
				fin = leaves[j]->fin;
				for (int i = debut; i < fin; i++) {
					Triangle t(vertices[indices[i].vtxi], vertices[indices[i].vtxj], vertices[indices[i].vtxk], Vector(1, 1, 1), false, false);
					Vector locP, locN, loccolor;
					if (t.intersect(r, locP, locN, loccolor)) {
						double dist = (locP - r.C).norm2();
						found = true;
						if (dist < dist_min) {
							dist_min = dist;
							P = locP;
							double alpha, beta, gamma;
							t.get_barycentric(alpha, beta, gamma, P);
							// N = alpha * normals[indices[i].ni] + beta * normals[indices[i].nj] + gamma * normals[indices[i].nk];
							N = locN;
							N.normalize();
							Vector uv = uvs[indices[i].uvi] * alpha + beta * uvs[indices[i].uvj]  + gamma * uvs[indices[i].uvk];
							int width = textures_width[indices[i].group];
							int height = textures_width[indices[i].group];
							int x = fabs(fmod(uv.x, 1.)) * width;
							int y = fabs(fmod(uv.y, 1.)) * height;
							color.x = textures[indices[i].group][(x + y * width) * 3] / 255.;
							color.y = textures[indices[i].group][(x + y * width) * 3 + 1] / 255.;
							color.z = textures[indices[i].group][(x + y * width) * 3 + 2] / 255.;
						}

					}
				}
			}
			return found;
		}
	}
	void add_textures(const char* filename) {
		int w, h, c;
		textures.push_back(stbi_load(filename, &w, &h, &c, 3));
		textures_width.push_back(w);
		textures_height.push_back(h);
	}
	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	std::vector<unsigned char*> textures;
	std::vector<int> textures_width;
	std::vector<int> textures_height;
	BBox bb;
	Vector offset;
	double sampling;
	BVHNode geometry_node;
};

class Scene {
public:
	Scene(int I, Sphere& L) : I(I), L(L) {};
	void addObject(Object* s) {
		objects.push_back(s);
	}
	bool intersect(const Ray& r, Vector& P, Vector& N, Object* &S, Vector& color) {
		double dist_min = 1E15;
		bool found = false;
		for (int i = 0; i < objects.size(); i++) {
			Object* Sloc = objects[i];
			Vector Ploc, Nloc, colorloc;
			if (Sloc->intersect(r, Ploc, Nloc, colorloc)) {
				double dist_cam = (Ploc - r.C).norm2();
				if (dist_cam < dist_min) {
					dist_min = dist_cam;
					P = Ploc;
					N = Nloc;
					S = Sloc;
					color = colorloc;
				}
				found = true;
			}
		}
		return found;
	}

	Vector create_T1(const Vector& N)
	{
		Vector T1;
		if (abs(N.x) <= abs(N.y) && abs(N.x) <= abs(N.z))
		{
			T1 = Vector(0., -N.z, N.y);
		}
		else
		{
			if (abs(N.y) <= abs(N.x) && abs(N.y) <= abs(N.z))
			{
				T1 = Vector(-N.z, 0., N.x);
			}
			else
			{
				T1 = Vector(-N.y, N.x, 0.);
			}
		}
		T1.normalize();
		return T1;
	}

	Vector getRandomPointLight(Sphere& L, Vector& P)
	{
		Vector PL = P - L.O;
		PL.normalize();
		double r1 = U(e);
		double r2 = U(e);
		Vector T1 = create_T1(PL);
		Vector T2 = cross(PL, T1);
		Vector wi_1 = cos(2 * M_PI*r1)*sqrt(1 - r2) * T1;
		Vector wi_2 = sin(2 * M_PI*r1)*sqrt(1 - r2) * T2;
		Vector wi_3 = sqrt(r2) * PL;
		Vector direction = wi_1 + wi_2 + wi_3;
		direction.normalize();
		Vector xPrime = L.O + direction * L.R;
		return xPrime;
	}

	Vector getColor(Ray& r, double& epsilon, int bounce, double n_sphere, double n_air)
	{
		Vector pixelColor(0, 0, 0);
		Vector P, N, color;
		Object* S;
		if (intersect(r, P, N, S, color) && bounce > 0) {
			if (S->mirror)
			{
				bounce = bounce - 1;
				Ray r_reflect(P + epsilon * N, r.u - 2.*dot(r.u, N)*N);
				pixelColor = getColor(r_reflect, epsilon, bounce, n_sphere, n_air);
			}
			else if (S->transparent)
			{
				double dotP = dot(r.u, N);
				Vector r_ref_v;
				if (dotP > 0) {
					Vector r_ref_t = n_sphere / n_air * (r.u - dotP * N);
					Vector r_ref_v = r_ref_t + sqrt(1 - pow(n_sphere / n_air, 2.)*(1 - pow(dotP, 2.)))*N;
					Ray r_ref = Ray(P + epsilon * N, r_ref_v);
					pixelColor = getColor(r_ref, epsilon, bounce, n_sphere, n_air);
				}
				else 
				{
					Vector r_ref_t = n_air / n_sphere * (r.u - dotP * N);
					Vector r_ref_v = r_ref_t - sqrt(1 - pow(n_air / n_sphere, 2.)*(1 - pow(dotP, 2.)))*N;
					Ray r_ref = Ray(P - epsilon * N, r_ref_v);
					pixelColor = getColor(r_ref, epsilon, bounce, n_sphere, n_air);
				}
			}
			else if (S->light)
			{
				pixelColor = Vector(0, 0, 0);
			}
			else
			{
				bounce = bounce - 1;
				Vector xPrime = getRandomPointLight(L, P);
				Vector PxPrime = xPrime - P;
				double dist_lum = PxPrime.norm2();
				PxPrime.normalize();
				Vector Nprimee = xPrime - L.O;
				Nprimee.normalize();
				Vector OX = P - L.O;
				OX.normalize();
				Vector Prime, Nprime, colorprime;
				Object* Sprime;
				Vector Pespsilon = P + epsilon * PxPrime; // to avoid noisy image
				Ray rprime = Ray(Pespsilon, PxPrime);
				if (intersect(rprime, Prime, Nprime, Sprime, colorprime)) { // shadow detection
					Vector PrimeP = P - Prime;
					if (PrimeP.norm2() < dist_lum)
					{
						pixelColor = Vector(0, 0, 0);
					}
					else
					{
						pixelColor = I / (4 * M_PI*dist_lum) * (color / M_PI) * 
							dot(PxPrime, N) * dot(-1. * PxPrime, Nprimee) / dot(Nprimee, OX);
					}
				}
				else
				{
					pixelColor = I / (4 * M_PI*dist_lum) * (color / M_PI) * dot(PxPrime, N) * 
						dot(-1. * PxPrime, Nprimee) / dot(Nprimee, OX);
				}
				double r1 = U(e);
				double r2 = U(e);
				Vector T1 = create_T1(N);
				Vector T2 = cross(N, T1);
				Vector wi_1 = cos(2 * M_PI*r1)*sqrt(1 - r2) * T1;
				Vector wi_2 = sin(2 * M_PI*r1)*sqrt(1 - r2) * T2;
				Vector wi_3 = sqrt(r2) * N;
				Ray r_diffusion = Ray(P + epsilon * N, wi_1 + wi_2 + wi_3);
				pixelColor = pixelColor + color * getColor(r_diffusion, epsilon, bounce, n_sphere, n_air);
			}
		}
		return pixelColor;
	};
	std::vector<Object*> objects;
	int I;
	Sphere L;
};

int main() {
	int W = 512;
	int H = 512;

	// Camera
	double fov = 60 * M_PI / 180.;
	Vector C(0, 0, 55);
	std::vector<unsigned char> image(W*H * 3, 0);
	double focal_distance = 55.; 

	// Light
	Sphere L(Vector(-10, 20, 40), 10, Vector(1., 1., 1.), false, false);
	L.light = true;
	double I = 214000000;

	// Noise
	double epsilon = 0.01;

	// Sphere indice
	double n_sphere = 1.4;
	double n_air = 1.;

	Sphere s1(Vector(-15, 0, 0), 5, Vector(1., 1., 1.), true, false);
	Sphere s2(Vector(15, 0, 0), 5, Vector(0., 1., 0.), false, true);
	Sphere s3(Vector(0, -1000, 0), 1000 - 10, Vector(0., 0., 1.), false, false); //sol
	Sphere s4(Vector(0, 0, -1050), 1000 - 60, Vector(1., 0., 1.), false, false); //fond
	Sphere s5(Vector(0, 0, 1000), 1000 - 60, Vector(1., 1., 0.), false, false); //invisible
	Sphere s6(Vector(1000, 0, 0), 1000 - 60, Vector(1., 0., 0.), false, false); //droite
	Sphere s7(Vector(-1000, 0, 0), 1000 - 60, Vector(0., 1., 0.), false, false); //gauche
	Sphere s8(Vector(0, 1000, 0), 1000 - 60, Vector(1., 1., 1.), false, false); //plafond
	Geometry g1(Vector(0., -5., 5.), 15, Vector(1, 1, 1), false, false);
	g1.readOBJ("cadnav.com_model/Model_D0515019/Beautiful Girl.obj", false);
	g1.add_textures("cadnav.com_model/Model_D0515019/visage.bmp");
	g1.add_textures("cadnav.com_model/Model_D0515019/cheveux.bmp");
	g1.add_textures("cadnav.com_model/Model_D0515019/corps.bmp");
	g1.add_textures("cadnav.com_model/Model_D0515019/pantalon.bmp");
	g1.add_textures("cadnav.com_model/Model_D0515019/accessoires.bmp");
	g1.add_textures("cadnav.com_model/Model_D0515019/mains.bmp");
	Scene scene(I, L);
	scene.addObject(&s1);
	scene.addObject(&s2);
	scene.addObject(&g1);
	scene.addObject(&s3);
	scene.addObject(&s4);
	scene.addObject(&s5);
	scene.addObject(&s6);
	scene.addObject(&s7);
	scene.addObject(&s8);
	scene.addObject(&L);



	// Creating the pixels
	for (int i = 0; i < W; i++) {
		cout << i << endl;
		for (int j = 0; j < H; j++) {
			Vector pixelColor;
			for (int k=0; k<30; k++)
			{
				double r1 = U(e);
				double r2 = U(e);
				double r3 = U(e);
				double r4 = U(e);
				double offsetx = cos(2 * M_PI*r1)*sqrt(-2*log(r2))*0.5;
				double offsety = sin(2 * M_PI*r1)*sqrt(-2 * log(r2))*0.5;
				double offsetx_C = cos(2 * M_PI*r3)*sqrt(-2 * log(r4));
				double offsety_C = sin(2 * M_PI*r3)*sqrt(-2 * log(r4));
				Vector CPrime = Vector(C.x + offsetx_C, C.y + offsety_C, C.z);
				Vector u = Vector(j - W / 2 + offsetx, -i + H / 2 + offsety, -W / (2 * tan(fov / 2)));
				u.normalize();
				Vector focal_P = C + focal_distance * u;
				Vector uPrime = focal_P - CPrime;
				uPrime.normalize();
				Ray r = Ray(CPrime, uPrime);
				pixelColor = pixelColor + scene.getColor(r, epsilon, 5, n_sphere, n_air);
			}
			image[(i*W + j) * 3 + 0] = std::min(255., std::pow(pixelColor.x, 0.45));
			image[(i*W + j) * 3 + 1] = std::min(255., std::pow(pixelColor.y, 0.45));
			image[(i*W + j) * 3 + 2] = std::min(255., std::pow(pixelColor.z, 0.45));
		}
	}
	cout << "Writing image" << endl;
	stbi_write_png("image.png", W, H, 3, &image[0], 0);
	cout << "Finished writing image" << endl;
	return 0;
}
