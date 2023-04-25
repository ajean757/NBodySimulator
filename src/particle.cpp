#include <nanogui/nanogui.h>

#include "misc/sphere_drawing.h"
#include "particle.h"

using namespace nanogui;
using namespace CGL;

void Particle::render(GLShader& shader) {
	m_sphere_mesh.draw_sphere(shader, origin, radius);
}