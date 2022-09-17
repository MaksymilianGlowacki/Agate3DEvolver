#ifndef AGATE3DEVOLVER_MESH_FUNCTIONS_FMM_H
#define AGATE3DEVOLVER_MESH_FUNCTIONS_FMM_H

#include "../mesh_functions.h"
#include <vector>
#include <algorithm>
#include <array>

std::vector<std::array<int32_t, 3>> find_seeds(triangle *mesh, unsigned long X, unsigned long Y, unsigned long Z,  size_t mesh_size);
void find_triangle_grid_points(bool *grid_points, triangle t, unsigned long X, unsigned long Y, unsigned long Z);
bool within_range(triangle t, double px, double py, double pz);
bool outside_mesh(triangle *mesh, int px, int py, int pz, size_t mesh_size, int Y);
std::vector<std::array<size_t, 3>> find_seeds_fim(triangle *mesh, unsigned long X, unsigned long Y, unsigned long Z,  size_t mesh_size);



#endif //AGATE3DEVOLVER_MESH_FUNCTIONS_FMM_H
