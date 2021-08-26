
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cfloat> // DBL_MAX

#ifdef __APPLE__
/* Defined before OpenGL and GLUT includes to avoid deprecation messages */
#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define epsilon 1.0e-5

struct Ray {
	double o[3], d[3]; // origin, direction
};

struct Colourtype {
	float r, g, b;
};

void crossProductd(double* cp, double* a, double* b) {
	cp[0] = a[1] * b[2] - a[2] * b[1];
	cp[1] = a[2] * b[0] - a[0] * b[2];
	cp[2] = a[0] * b[1] - a[1] * b[0];
}

void addcolour(Colourtype& Colour, Colourtype& toadd) {
	Colour.r += toadd.r;
	Colour.g += toadd.g;
	Colour.b += toadd.b;
}

void multColour(Colourtype& Colour, float scalar) {
	Colour.r *= scalar;
	Colour.g *= scalar;
	Colour.b *= scalar;
}

void clampColour(Colourtype& Colour) {
	if (Colour.r > 1.0f) Colour.r = 1.0f;
	if (Colour.g > 1.0f) Colour.g = 1.0f;
	if (Colour.b > 1.0f) Colour.b = 1.0f;
}

double normalize3d(double* v) {
	double f = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if (f > 0.0) {
		f = sqrt(f);
		double f1 = 1.0 / f;
		v[0] *= f1; v[1] *= f1; v[2] *= f1;
	}
	return f;
}


struct cameratype {
	double eye[3], dir[3], up[3];
	double zNear, fovy;
};

enum objtype { sphere = 0, plane, triangle, cylinder, torus };

struct Objecttype {
	objtype type;
	Colourtype Ra, Rd, Rs; // Ambient, diffuse, and specular coefficients
	float f; // specular highlight coefficient
	double n, rc, tc; // index of refraction, reflection coefficient, transmission coefficient
	bool has_refl, has_trans;
	double pdata[12]; // points, normals, different for each objtype
	double getIntersection(Ray* rvec, double* pt, double* norm);
	double sphere_intersect(Ray* ray, double* pt, double* norm);
	double plane_intersect(Ray* ray, double* pt);
	double triangle_intersect(Ray* ray, double* pt);
};


double Objecttype::getIntersection(Ray* rvec, double* pt, double* norm) {
	switch (type) {
	case sphere:
		return sphere_intersect(rvec, pt, norm);
		break;
	case plane:
		norm[0] = pdata[0]; norm[1] = pdata[1]; norm[2] = pdata[2];
		return plane_intersect(rvec, pt);
		break;
	case triangle:
		norm[0] = pdata[9]; norm[1] = pdata[10]; norm[2] = pdata[11];
		return triangle_intersect(rvec, pt);
		break;
	}
}

double Objecttype::sphere_intersect(Ray* ray, double* pt, double* norm) {
	double vx = (pdata[0] - ray->o[0]);
	double vy = (pdata[1] - ray->o[1]);
	double vz = (pdata[2] - ray->o[2]);
	double d = vx * ray->d[0] + vy * ray->d[1] + vz * ray->d[2];
	double t = pdata[3] * pdata[3] + d * d - vx * vx - vy * vy - vz * vz;
	if (t < epsilon * epsilon) return DBL_MAX; //no intersection
	t = sqrt(t);
	if (d - t < epsilon) {
		if (d < 0.0) return DBL_MAX;
		t += d;
	}
	else {
		t = d - t;
	}
	if (pt) {
		d = 1.0 / pdata[3];
		for (int i = 0; i != 3; i++) {
			pt[i] = t * ray->d[i] + ray->o[i];
			norm[i] = d * (pt[i] - pdata[i]);
		}
	}
	return t;
}

double Objecttype::plane_intersect(Ray* ray, double* pt) {
	double d = pdata[0] * ray->d[0] + pdata[1] * ray->d[1] + pdata[2] * ray->d[2];
	if (fabs(d) < epsilon)
		return DBL_MAX;
	double t = -(pdata[0] * ray->o[0] + pdata[1] * ray->o[1] + pdata[2] * ray->o[2] + pdata[3]) / d;
	if (t > epsilon) {
		if (pt)
			for (int i = 0; i != 3; i++)
				pt[i] = ray->o[i] + t * ray->d[i];

		return t;
	}
	return DBL_MAX;
}

// modified from https://en.m.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
double Objecttype::triangle_intersect(Ray* ray, double* pt) {
	//	Vector3D edge1, edge2, h, s, q;
	double t, u, v;
	double h[3], s[3], q[3];
	//	edge1 = vertex1 - vertex0;
	//	edge2 = vertex2 - vertex0;
	crossProductd(h, ray->d, pdata + 6);
	//	h = rayVector.crossProduct(edge2);
	t = pdata[3] * h[0] + pdata[4] * h[1] + pdata[5] * h[2];
	//	a = edge1.dotProduct(h);
	if (fabs(t) < epsilon)
		return DBL_MAX;	// This ray is parallel to this triangle.
	t = 1.0 / t;
	for (int i = 0; i != 3; i++)
		s[i] = ray->o[i] - pdata[i];
	u = t * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
	//	s = rayOrigin - vertex0;
	//	u = f * s.dotProduct(h);
	if (u < 0.0 || u > 1.0)
		return DBL_MAX;
	crossProductd(q, s, pdata + 3);
	//	q = s.crossProduct(edge1);
	v = t * (ray->d[0] * q[0] + ray->d[1] * q[1] + ray->d[2] * q[2]);
	//	v = f * rayVector.dotProduct(q);
	if (v < 0.0 || u + v > 1.0)
		return DBL_MAX;
	// At this stage we can compute t to find out where the intersection point is on the line.
	t *= pdata[6] * q[0] + pdata[7] * q[1] + pdata[8] * q[2];
	//	float t = f * edge2.dotProduct(q);
	if (t > epsilon) // ray intersection
	{
		if (pt)
			for (int i = 0; i != 3; i++)
				pt[i] = ray->o[i] + t * ray->d[i];
		//outIntersectionPoint = rayOrigin + rayVector * t;
		return t;
	}
	// This means that there is a line intersection but not a ray intersection.
	return DBL_MAX;
}

struct Lighttype {
	double pos[3]; // light position
	Colourtype Ia, Is; // Ambient light intensity and light intensity
	//bool active;
};

// Global variables
Colourtype background_colour;
std::vector <Lighttype> Light;
std::vector <Objecttype> Object;
cameratype cam;
std::vector <char> imagedata;
int imagewidth, imageheight;
int MAX_DEPTH;

bool RayIntersection(Ray* rvec, Objecttype** obj, double* pt, double* norm) {
	double mind = DBL_MAX, d, p[3], n[3];
	*obj = 0;
	for (size_t i = 0; i != Object.size(); i++) {
		d = Object[i].getIntersection(rvec, p, n);
		if (d < mind) {
			mind = d;
			pt[0] = p[0]; pt[1] = p[1]; pt[2] = p[2];
			norm[0] = n[0]; norm[1] = n[1]; norm[2] = n[2];
			*obj = &Object[i];
		}
	}
	return mind != DBL_MAX;
}

Colourtype shade(double* R, Objecttype* obj, double* pt, double* norm) {
	Colourtype Colour = { 0.0f, 0.0f, 0.0f };
	Colourtype c;
	for (size_t i = 0; i != Light.size(); i++) {
		Ray lray;
		c.r = Light[i].Ia.r * obj->Ra.r;
		c.g = Light[i].Ia.g * obj->Ra.g;
		c.b = Light[i].Ia.b * obj->Ra.b;

		for (int j = 0; j != 3; j++) {
			lray.o[j] = pt[j];
			lray.d[j] = Light[i].pos[j] - pt[j];
		}
		double d = normalize3d(lray.d);
		bool lighting = true;
		// test if there are other objects between the light and the object.
		for (size_t j = 0; j != Object.size(); j++) {
			switch (Object[j].type) {
			case sphere:
				lighting = Object[j].sphere_intersect(&lray, 0, 0) > d;
				break;
			case plane:
				lighting = Object[j].plane_intersect(&lray, 0) > d;
				break;
			case triangle:
				lighting = Object[j].triangle_intersect(&lray, 0) > d;
				break;
			}
			if (!lighting)
				break;
		}
		if (lighting) {
			double cosA = norm[0] * lray.d[0] + norm[1] * lray.d[1] + norm[2] * lray.d[2];
			if (cosA > 0.0) {
				c.r += cosA * (Light[i].Is.r * obj->Rd.r);
				c.g += cosA * (Light[i].Is.g * obj->Rd.g);
				c.b += cosA * (Light[i].Is.b * obj->Rd.b); // diffuse contribution
				// Combining reflection and specular calculations (CosA now is cos(B)):
				cosA = R[0] * lray.d[0] + R[1] * lray.d[1] + R[2] * lray.d[2];
				if (cosA > 0.0) { // specular contribution
					cosA = pow(cosA, obj->f);
					c.r += cosA * (Light[i].Is.r * obj->Rs.r);
					c.g += cosA * (Light[i].Is.g * obj->Rs.g);
					c.b += cosA * (Light[i].Is.b * obj->Rs.b);
				}
			}
		}
		addcolour(Colour, c);
	}
	return Colour;
}

/*void calc_transmission(Ray *rtr, Ray *rvec, Objecttype *obj, double *pt, double *norm) {
}*/

Colourtype RayTrace(Ray* rvec, int rec_depth) {
	if (rec_depth > MAX_DEPTH)
		return background_colour;
	double hitpt[3], hitnorm[3];
	Objecttype* hitobj;
	if (RayIntersection(rvec, &hitobj, hitpt, hitnorm)) {
		Ray t;
		double Ci = -hitnorm[0] * rvec->d[0] - hitnorm[1] * rvec->d[1] - hitnorm[2] * rvec->d[2];
		double R[3], d = 2.0 * Ci, w = 1.0;
		for (int i = 0; i != 3; i++)
			R[i] = rvec->d[i] + d * hitnorm[i]; //  Used in reflection and specular calculations
		Colourtype Colour, local_colour = shade(R, hitobj, hitpt, hitnorm);
		if (hitobj->has_refl) {
			//calc_reflection(&t, rvec, hitobj, hitpt, hitnorm);
			for (int i = 0; i != 3; i++) {
				t.d[i] = R[i];
				t.o[i] = hitpt[i];
			}
			w -= hitobj->rc;
			Colour = RayTrace(&t, rec_depth + 1);
			multColour(Colour, hitobj->rc);
		}
		else
			Colour.r = Colour.g = Colour.b = 0.0f;
		if (hitobj->has_trans) {
			//calc_transmission(&t, rvec, hitobj, hitpt, hitnorm);
			double n = (Ci < 0 ? hitobj->n : 1.0 / hitobj->n);
			double b = (1.0 + n * n * (Ci * Ci - 1.0));
			if (b >= 0.0) {
				b = n * Ci + (Ci < 0 ? sqrt(b) : -sqrt(b));
				//b = n*Ci - sqrt(b);
				for (int i = 0; i != 3; i++) {
					t.d[i] = n * rvec->d[i] + b * hitnorm[i];
					t.o[i] = hitpt[i];
				}
				w -= hitobj->tc;
				Colourtype c = RayTrace(&t, rec_depth + 1);
				multColour(c, hitobj->tc);
				addcolour(Colour, c);
			}
		}
		multColour(local_colour, w);
		addcolour(Colour, local_colour);
		return Colour;
	}
	return background_colour;
}

void start() {
	imagedata.resize(3 * imagewidth * imageheight);
	double ux[3], uy[3], o[3], ratio, t, dx, dy;
	Ray rvec;
	Colourtype Colour;
	int k = -1;
	normalize3d(cam.dir);
	crossProductd(ux, cam.dir, cam.up);
	normalize3d(ux);
	crossProductd(uy, ux, cam.dir);
	t = cam.zNear * tan(cam.fovy * (atan(1.0) / 45));
	ratio = (double)imagewidth / imageheight;
	for (int i = 0; i != 3; i++)
		o[i] = cam.eye[i] + cam.zNear * cam.dir[i] - t * (ratio * ux[i] + uy[i]);
	dx = ratio * 2 * t / (imagewidth - 1);
	dy = 2 * t / (imageheight - 1);
	for (int i = 0; i != 3; i++) {
		ux[i] *= dx;
		uy[i] *= dy;
	}
	for (int j = 0; j != imageheight; j++)
		for (int i = 0; i != imagewidth; i++) {
			for (int h = 0; h != 3; h++) {
				rvec.o[h] = o[h] + i * ux[h] + j * uy[h];
				rvec.d[h] = rvec.o[h] - cam.eye[h];
			}
			normalize3d(rvec.d);
			Colour = RayTrace(&rvec, 0);
			clampColour(Colour);
			imagedata[++k] = 255 * Colour.r;
			imagedata[++k] = 255 * Colour.g;
			imagedata[++k] = 255 * Colour.b;
		}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(imagewidth, imageheight, GL_RGB, GL_UNSIGNED_BYTE, &imagedata[0]);
	glFlush();
}

void fileerror(int line) {
	std::cout << "Invalid value on line " << line << "\n";
	imagewidth = -1;
}

void readSceneFile(const char* filename) {
	imagewidth = -1;
	std::ifstream in(filename);
	if (!in)
		return;
	Lighttype lgt;
	Objecttype obj;
	char cs[256];
	int line = 0;
	MAX_DEPTH = 5;
	background_colour.r = background_colour.g = background_colour.b = 0.0f;
	// default object properties:
	obj.has_refl = obj.has_trans = false;
	obj.Rs.r = obj.Rs.g = obj.Rs.b = 0.9f;
	obj.Ra.r = obj.Rd.r = obj.Ra.g = obj.Rd.g = obj.Ra.b = obj.Rd.b = 0.4f;
	obj.f = 10.0;
	obj.n = 1.15;
	// default camera settings:
	cam.up[0] = cam.up[2] = cam.dir[0] = cam.dir[1] = cam.dir[2] = cam.eye[0] = cam.eye[1] = 0.0;
	cam.eye[2] = 3.0;
	cam.up[1] = 1.0;
	cam.zNear = 0.1;
	cam.fovy = 45.0;
	// default light properties:
	lgt.pos[0] = lgt.pos[1] = lgt.pos[2] = 1.0;
	lgt.Ia.r = lgt.Ia.g = lgt.Ia.b = 0.1f;
	lgt.Is.r = lgt.Is.g = lgt.Is.b = 1.0f;
	in.getline(cs, 255);
	while (in.good()) {
		line++;
		//std::cout << cs << "\n";
		if (!strncmp(cs, "width:", 6)) {
			imagewidth = atoi(cs + 6);
			if (imagewidth % 4)
				imagewidth = 4 * ((imagewidth + 3) / 4); // to avoid problem with glDrawPixels and the bmp file format
			if (imagewidth < 1) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "height:", 7)) {
			imageheight = atoi(cs + 7);
			if (imageheight < 1) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "max recursion level:", 20)) {
			MAX_DEPTH = atoi(cs + 20);
			if (MAX_DEPTH < 1) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "background colour:", 18)) {
			char* c;
			background_colour.r = strtod(cs + 18, &c);
			background_colour.g = strtod(c, &c);
			background_colour.b = strtod(c, 0);
		}
		else if (!strncmp(cs, "z near:", 7)) {
			cam.zNear = atof(cs + 7);
			if (cam.zNear <= 0.0) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "fovy:", 5)) {
			cam.fovy = atof(cs + 5);
			if (cam.fovy <= 0.0 || cam.fovy > 270.0) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "camera position:", 16)) {
			char* c;
			cam.eye[0] = strtod(cs + 16, &c);
			cam.eye[1] = strtod(c, &c);
			cam.eye[2] = strtod(c, 0);
			if (*c == 0) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "target position:", 16)) {
			char* c;
			cam.dir[0] = strtod(cs + 16, &c);
			cam.dir[1] = strtod(c, &c);
			cam.dir[2] = strtod(c, 0);
			if (*c == 0) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "up vector:", 10)) {
			char* c;
			cam.up[0] = strtod(cs + 10, &c);
			cam.up[1] = strtod(c, &c);
			cam.up[2] = strtod(c, 0);
			if (*c == 0) {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "light:", 6)) {
			char* c;
			lgt.pos[0] = strtod(cs + 6, &c);
			lgt.pos[1] = strtod(c, &c);
			lgt.pos[2] = strtod(c, &c);
			if (strlen(c) > 5) {
				lgt.Ia.r = strtod(c, &c);
				lgt.Ia.g = strtod(c, &c);
				lgt.Ia.b = strtod(c, &c);
				if (strlen(c) > 5) {
					lgt.Is.r = strtod(c, &c);
					lgt.Is.g = strtod(c, &c);
					lgt.Is.b = strtod(c, 0);
				}
			}
			Light.push_back(lgt);
		}
		else if (!strncmp(cs, "object", 6)) {
			char* c = cs + 6;
			if (!strncmp(c, " property", 9)) {
				c += 9;
				if (!strncmp(c, " colour:", 8)) {
					obj.Ra.r = strtod(c + 8, &c);
					obj.Ra.g = strtod(c, &c);
					obj.Ra.b = strtod(c, &c);
					if (strlen(c) > 5) {
						obj.Rd.r = strtod(c, &c);
						obj.Rd.g = strtod(c, &c);
						obj.Rd.b = strtod(c, &c);
						if (strlen(c) > 5) {
							obj.Rs.r = strtod(c, &c);
							obj.Rs.g = strtod(c, &c);
							obj.Rs.b = strtod(c, 0);
						}
						else
							obj.Rs.r = obj.Rs.g = obj.Rs.b = 0.0f;
					}
					else {
						obj.Rd.r = obj.Ra.r;
						obj.Rd.g = obj.Ra.g;
						obj.Rd.b = obj.Ra.b;
						obj.Rs.r = obj.Rs.g = obj.Rs.b = 0.0f;
					}
				}
				else if (!strncmp(c, " transmission:", 14)) {
					obj.tc = atof(c + 14);
					if (obj.tc < 0.0 || obj.tc > 1.0) {
						fileerror(line);
						break;
					}
					if (obj.tc + obj.rc > 1.0)
						obj.rc = 1.0 - obj.tc;
					obj.has_trans = obj.tc > 0.0;
				}
				else if (!strncmp(c, " reflection:", 12)) {
					obj.rc = atof(c + 12);
					if (obj.rc < 0.0 || obj.rc > 1.0) {
						fileerror(line);
						break;
					}
					if (obj.tc + obj.rc > 1.0)
						obj.tc = 1.0 - obj.rc;
					obj.has_refl = obj.rc > 0.0;
				}
				else if (!strncmp(c, " shininess:", 11)) {
					obj.f = atof(c + 11);
					if (obj.f < 1.0) {
						fileerror(line);
						break;
					}
				}
				else if (!strncmp(c, " refraction index:", 18)) {
					obj.n = atof(c + 18);
					if (obj.n < 1.0) {
						fileerror(line);
						break;
					}
				}
				else {
					fileerror(line);
					break;
				}
			}
			else if (!strncmp(c, " sphere:", 8)) {
				obj.type = sphere;
				obj.pdata[0] = strtod(c + 8, &c);
				obj.pdata[1] = strtod(c, &c);
				obj.pdata[2] = strtod(c, &c);
				obj.pdata[3] = strtod(c, 0);
				Object.push_back(obj);
				if (obj.pdata[3] <= 0.0) {
					fileerror(line);
					break;
				}
			}
			else if (!strncmp(c, " plane:", 7)) {
				obj.type = plane;
				obj.pdata[0] = strtod(c + 7, &c);
				obj.pdata[1] = strtod(c, &c);
				obj.pdata[2] = strtod(c, &c);
				if (normalize3d(obj.pdata) == 0.0) {
					fileerror(line);
					break;
				}
				obj.pdata[4] = strtod(c, &c);
				obj.pdata[5] = strtod(c, &c);
				obj.pdata[6] = strtod(c, 0);
				obj.pdata[3] = obj.pdata[0] * obj.pdata[4] + obj.pdata[1] * obj.pdata[5] + obj.pdata[2] * obj.pdata[6];
				Object.push_back(obj);
			}
			else if (!strncmp(c, " triangle:", 10)) {
				obj.type = triangle;
				obj.pdata[0] = strtod(c + 10, &c);
				obj.pdata[1] = strtod(c, &c);
				obj.pdata[2] = strtod(c, &c);
				obj.pdata[3] = strtod(c, &c);
				obj.pdata[4] = strtod(c, &c);
				obj.pdata[5] = strtod(c, &c);
				obj.pdata[6] = strtod(c, &c);
				obj.pdata[7] = strtod(c, &c);
				if (strlen(c) < 1) {
					fileerror(line);
					break;
				}
				obj.pdata[8] = strtod(c, 0);
				for (int i = 0; i != 3; i++) {
					obj.pdata[3 + i] -= obj.pdata[i]; // change point1 for edge1 vector
					obj.pdata[6 + i] -= obj.pdata[i]; // change point2 for edge2 vector
				}
				crossProductd(obj.pdata + 9, obj.pdata + 3, obj.pdata + 6); // normal
				normalize3d(obj.pdata + 9);
				Object.push_back(obj);
			}
			else {
				fileerror(line);
				break;
			}
		}
		else if (!strncmp(cs, "end", 3))
			break;
		else if (cs[0] != '#' && strlen(cs)) {
			fileerror(line);
			break;
		}
		in.getline(cs, 255);
	}
	for (int i = 0; i != 3; i++)
		cam.dir[i] -= cam.eye[i]; // now dir contains a vector from eye to target position
	if (Light.empty())
		Light.push_back(lgt);
	if (Object.empty()) {
		std::cout << "Objects not found!\n";
		imagewidth = -1;
	}
}
// modified from https://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
void savebmp(const char* filename, int w, int h, char* rgb) {
	std::string theFile = std::string(filename);
	std::string clean = theFile.substr(0, theFile.length() - 4); //removes weird .txt extension from bmp file
	std::string bmpname = std::string(clean) + ".bmp";
	std::ofstream stream(bmpname.c_str(), std::ios::binary);
	if (!stream)
		return;

	// mimeType = "image/bmp";
	unsigned char file[14] = {
		'B','M', // magic
		0,0,0,0, // size in bytes
		0,0, // app data
		0,0, // app data
		40 + 14,0,0,0 // start of data offset
	};
	unsigned char info[40] = {
		40,0,0,0, // info hd size
		0,0,0,0, // width
		0,0,0,0, // heigth
		1,0, // number color planes
		24,0, // bits per pixel
		0,0,0,0, // compression is none
		0,0,0,0, // image bits size
		0x13,0x0B,0,0, // horz resoluition in pixel / m
		0x13,0x0B,0,0, // vert resolutions (0x03C3 = 96 dpi, 0x0B13 = 72 dpi)
		0,0,0,0, // #colors in pallete
		0,0,0,0, // #important colors
	};

	int padSize = (4 - (w * 3) % 4) % 4;
	int sizeData = w * h * 3 + h * padSize;
	int sizeAll = sizeData + sizeof(file) + sizeof(info);

	file[2] = (unsigned char)(sizeAll);
	file[3] = (unsigned char)(sizeAll >> 8);
	file[4] = (unsigned char)(sizeAll >> 16);
	file[5] = (unsigned char)(sizeAll >> 24);

	info[4] = (unsigned char)(w);
	info[5] = (unsigned char)(w >> 8);
	info[6] = (unsigned char)(w >> 16);
	info[7] = (unsigned char)(w >> 24);

	info[8] = (unsigned char)(h);
	info[9] = (unsigned char)(h >> 8);
	info[10] = (unsigned char)(h >> 16);
	info[11] = (unsigned char)(h >> 24);

	info[20] = (unsigned char)(sizeData);
	info[21] = (unsigned char)(sizeData >> 8);
	info[22] = (unsigned char)(sizeData >> 16);
	info[23] = (unsigned char)(sizeData >> 24);

	stream.write((char*)file, sizeof(file));
	stream.write((char*)info, sizeof(info));

	char pad[3] = { 0,0,0 };

	int k = -1;
	for (int y = 0; y < h; y++)
	{
		for (int x = 0; x < w; x++)
		{
			char pixel[3];
			pixel[2] = rgb[++k];
			pixel[1] = rgb[++k];
			pixel[0] = rgb[++k];

			stream.write(pixel, 3);
		}
		stream.write(pad, padSize);
	}
}

// for test only
/*
void init() {
	imagewidth = 600;
	imageheight = 500;
	background_colour.r = background_colour.g = background_colour.b = 0.0f;

	cam.eye[0] = 0.0;
	cam.eye[1] = 3.0;
	cam.eye[2] = 7.0;
	cam.up[0] = cam.up[2] = 0.0;
	cam.up[1] = 1.0;
	cam.dir[0] = 0.0 - cam.eye[0];
	cam.dir[1] = 1.0 - cam.eye[1];
	cam.dir[2] = 0.0 - cam.eye[2];
	cam.zNear = 2;
	cam.fovy = 45;

	Lighttype lgt;
	lgt.pos[0] = -8;
	lgt.pos[1] = 20;
	lgt.pos[2] = 6;
	lgt.Ia.r = lgt.Ia.g = lgt.Ia.b = 0.1f;
	lgt.Is.r = lgt.Is.g = lgt.Is.b = 0.9f;
	Light.push_back(lgt);
	lgt.pos[0] = 8;
	lgt.pos[1] = 5;
	lgt.pos[2] = -4;
	lgt.Ia.r = lgt.Ia.g = lgt.Ia.b = 0.1f;
	lgt.Is.r = 0.9f; lgt.Is.g = lgt.Is.b = 0.6f;
	Light.push_back(lgt);

	Objecttype obj;
	obj.type = triangle;
	obj.Rs.r = obj.Rs.g = obj.Rs.b = 0.9f;
	obj.f = 20.0;
	obj.rc = 0.4;
	obj.has_refl = true;
	obj.has_trans = false;
	for (int j = -5; j != 5; j++)
		for (int i = -5; i != 5; i++) {
			obj.pdata[0] = i;
			obj.pdata[1] = 0;
			obj.pdata[2] = j;
			obj.pdata[3] = 0;
			obj.pdata[4] = 0;
			obj.pdata[5] = 1;
			obj.pdata[6] = 1;
			obj.pdata[7] = 0;
			obj.pdata[8] = 0;
			obj.pdata[9] = 0;
			obj.pdata[10] = 1;
			obj.pdata[11] = 0;
			if ((i + j)&1)
				obj.Ra.r = obj.Ra.g = obj.Ra.b = obj.Rd.r = obj.Rd.g = obj.Rd.b = 0.65f;
			else {
				obj.Ra.r = obj.Rd.r = 0.5f;
				obj.Ra.g = obj.Rd.g = obj.Ra.b = obj.Rd.b = 0.1f;
			}
			Object.push_back(obj);
			obj.pdata[0] = i + 1;
			obj.pdata[1] = 0;
			obj.pdata[2] = j + 1;
			obj.pdata[3] = 0;
			obj.pdata[4] = 0;
			obj.pdata[5] = -1;
			obj.pdata[6] = -1;
			obj.pdata[7] = 0;
			obj.pdata[8] = 0;
			Object.push_back(obj);
		}
	obj.type = sphere;
	obj.pdata[0] = 3;
	obj.pdata[1] = 2.5;
	obj.pdata[2] = -2.5;
	obj.pdata[3] = 2; // radius
	obj.rc = 0.1f;
	obj.Ra.r = obj.Rd.r = obj.Ra.b = obj.Rd.b = 0.1f;
	obj.Ra.g = obj.Rd.g = 0.6f;
	Object.push_back(obj);
	obj.pdata[0] = -3;
	obj.pdata[1] = 3;
	obj.pdata[2] = -2;
	obj.pdata[3] = 3; // radius
	obj.rc = 0.8f;
	obj.Ra.r = obj.Rd.r = obj.Ra.g = obj.Rd.g = 0.1f; obj.Ra.b = obj.Rd.b = 0.8f;
	Object.push_back(obj);
	obj.pdata[0] = 0;
	obj.pdata[1] = 2;
	obj.pdata[2] = 1;
	obj.pdata[3] = 2; // radius
	obj.Ra.r = obj.Ra.g = obj.Rd.r = obj.Rd.g = 0.7f;
	obj.Ra.b = obj.Rd.b = 0.9f;
	obj.has_trans = true;
	obj.n = 1.15;
	obj.rc = 0.2f;
	obj.tc = 0.7f;
	Object.push_back(obj);
}
*/

int main(int argc, char** argv)
{
	if (argc == 2)
		readSceneFile(argv[1]);
	else {
		std::cout << "Usage:\n\t" << argv[0] << "Scene_file\n";
		exit(1);
	}
	if (imagewidth == -1) {
		std::cout << "Error reading input file.\n";
		exit(1);
	}

	//	init();
	start();
	savebmp(argv[1], imagewidth, imageheight, &imagedata[0]);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(imagewidth, imageheight);
	glutInitWindowPosition(50, 50);
	glutCreateWindow("Ray Tracer");
	glutDisplayFunc(display);

	glutMainLoop();

	return 0;
}

