#ifndef PARTICLE_H
#define PARTICLE_H

#include "misc/sphere_drawing.h"
#include <CGL/vector3D.h>
#include "CGL/CGL.h"

using namespace CGL;
using namespace std;

struct Particle {
public:
	Particle(Vector3D origin, double radius, double mass) : origin(origin), radius(radius), mass(mass), m_sphere_mesh(Misc::SphereMesh()) {}

	void render(GLShader& shader);
	
private:
	Vector3D origin;
	double radius;
	double mass;
	Misc::SphereMesh m_sphere_mesh;
};

#endif /* PARTICLE_H */