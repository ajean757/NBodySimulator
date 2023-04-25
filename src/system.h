#ifndef SYSTEM_H
#define SYSTEM_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "particle.h"

using namespace CGL;
using namespace std;


struct System {
  System() {}
  ~System();

  void buildSystem();

  void simulate(double frames_per_sec, double simulation_steps, vector<Vector3D> external_accelerations/*,
    vector<CollisionObject*>* collision_objects*/);

  void reset();

  // void self_collide(PointMass& pm, double simulation_steps);

  // System properties
  //double width;
  //double height;
  //int num_width_points;
  //int num_height_points;
  //double thickness;
  //e_orientation orientation;

  // System components
  vector<Particle*> particles;
  // vector<vector<int>> pinned; might wanna keep this pinned list in case we wanna emulate light a solar system with a static sun in the middle or something?

};

#endif /* SYSTEM_H */
