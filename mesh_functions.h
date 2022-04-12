#ifndef AGATE3DEVOLVER_MESH_FUNCTIONS_H
#define AGATE3DEVOLVER_MESH_FUNCTIONS_H

#include <vector>
#include <cmath>

typedef long double triangle[3][3];

long double distance_to_triangle(triangle t, long double px, long double py, long double pz, size_t& inter);
long double distance_to_mesh(triangle *mesh, int x, int y, int z, size_t mesh_size);


#endif //AGATE3DEVOLVER_MESH_FUNCTIONS_H
