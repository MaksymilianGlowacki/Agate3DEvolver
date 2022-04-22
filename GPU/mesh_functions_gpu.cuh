#ifndef AGATE3DEVOLVER_MESH_FUNCTIONS_H
#define AGATE3DEVOLVER_MESH_FUNCTIONS_H

#include <vector>
#include <cmath>
#include <iostream>

typedef double triangle[3][3];

__device__
double distance_to_triangle(triangle t, double px, double py, double pz); //, size_t& inter);

__device__
double distance_to_mesh(triangle *mesh, double x, double y, double z, int mesh_size);


#endif //AGATE3DEVOLVER_MESH_FUNCTIONS_H
