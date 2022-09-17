#ifndef AGATE3DEVOLVER_MESH_FUNCTIONS_H
#define AGATE3DEVOLVER_MESH_FUNCTIONS_H

#include <vector>
#include <cmath>

typedef double triangle[3][3];

double distance_to_triangle(triangle t, double px, double py, double pz, size_t &inter, int Y);
double distance_to_mesh(triangle *mesh, int x, int y, int z, size_t mesh_size, int Y);


#endif //AGATE3DEVOLVER_MESH_FUNCTIONS_H
