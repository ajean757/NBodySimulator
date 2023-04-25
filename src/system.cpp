#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "system.h"

using namespace std;

System::~System() {
  // Destructor

  particles.clear();
}

void System::buildSystem() {
  Particle* p0 = new Particle(Vector3D(0.0), 1.0, 10.0, false);
  Particle* p1 = new Particle(Vector3D(0.5), 0.5, 10.0, false);
  particles = vector<Particle*>();
  particles.push_back(p0);
  particles.push_back(p1);
}

void System::simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations) {
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO: Compute net force on each particle and update positions accordingly
  particles[0]->position += Vector3D(0.001);

}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void System::reset() {
  Particle* p = particles[0];
  for (int i = 0; i < particles.size(); i++) {
    p->position = p->start_position;
    p->last_position = p->start_position;
    p++;
  }
}
