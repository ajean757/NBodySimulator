#ifndef PARTICLE_H
#define PARTICLE_H

#include "misc/sphere_drawing.h"
#include <CGL/vector3D.h>
#include "CGL/CGL.h"

using namespace CGL;
using namespace std;

struct Particle {
	Particle(Vector3D position, double radius, double mass, bool pinned)
		: position(position), radius(radius), mass(mass), pinned(pinned),
		start_position(position), m_sphere_mesh(Misc::SphereMesh()) {}
	
	void render(GLShader& shader);

	// static values
	double radius;
	double mass;
	bool pinned;
	Vector3D start_position;

	// dynamic values
	Vector3D position;
	Vector3D last_position;
	Vector3D forces;

private:
	Misc::SphereMesh m_sphere_mesh;

};

#endif /* PARTICLE_H */