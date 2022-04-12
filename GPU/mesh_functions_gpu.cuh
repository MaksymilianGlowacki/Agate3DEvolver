#ifndef AGATE3DEVOLVER_MESH_FUNCTIONS_H
#define AGATE3DEVOLVER_MESH_FUNCTIONS_H

#include <vector>
#include <cmath>

typedef double triangle[3][3];

__device__
double distance_to_triangle(triangle t, double px, double py, double pz, size_t& inter);

__device__
double distance_to_mesh(triangle *mesh, unsigned int x, unsigned int y, unsigned int z, size_t mesh_size);


#endif //AGATE3DEVOLVER_MESH_FUNCTIONS_H
