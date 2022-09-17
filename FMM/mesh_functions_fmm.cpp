#include "mesh_functions_fmm.h"

std::vector<std::array<int32_t, 3>> find_seeds(triangle *mesh, unsigned long X, unsigned long Y, unsigned long Z, size_t mesh_size) {
    auto grid_points = new bool[X * Y * Z];
    for(int i = 0; i < mesh_size; ++i)
    {
        find_triangle_grid_points(grid_points, mesh[i], X, Y, Z);
    }
    std::vector<std::array<int32_t, 3>> result;
    for(int z = 0; z < Z; ++z)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int x = 0; x < X; ++x)
            {
                if(grid_points[z * (X * Y) + y * X + x])
                {
                    result.emplace_back(std::array<int32_t, 3>({x, y, z}));
                }
            }
        }
    }
    delete[] grid_points;
    return result;
}

void find_triangle_grid_points(bool *grid_points, triangle t, unsigned long X, unsigned long Y, unsigned long Z)
{
    const int min_x = (int)std::floor(std::min({t[0][0], t[1][0], t[2][0]}));
    const int max_x = (int)std::ceil(std::max({t[0][0], t[1][0], t[2][0]}));

    const int min_y = (int)std::floor(std::min({t[0][1], t[1][1], t[2][1]}));
    const int max_y = (int)std::ceil(std::max({t[0][1], t[1][1], t[2][1]}));

    const int min_z = (int)std::floor(std::min({t[0][2], t[1][2], t[2][2]}));
    const int max_z = (int)std::ceil(std::max({t[0][2], t[1][2], t[2][2]}));
    for(int z = min_z; z <= max_z; ++z)
    {
        for(int y = min_y; y <= max_y; ++y)
        {
            for(int x = min_x; x <= max_x; ++x)
            {
                if((x >= 0) && (x < X) && (y >= 0) && (y < Y) && (z >= 0) && (z < Z))
                {
                    if(within_range(t, x, y, z))
                    {
                        grid_points[z * (X * Y) + y * X + x] = true;
                    }
                }
            }
        }
    }
}

bool within_range(triangle t, double px, double py, double pz)
{
    const double    x1 = t[0][0] - px, y1 = t[0][1] - py, z1 = t[0][2] - pz,
            x2 = t[1][0] - t[0][0], y2 = t[1][1] - t[0][1], z2 = t[1][2] - t[0][2],
            x3 = t[2][0] - t[0][0], y3 = t[2][1] - t[0][1], z3 = t[2][2] - t[0][2];
    double a, b, X = t[0][0], Y = t[0][1], Z = t[0][2];

    double d1 = (x1 * x1 +
                 y1 * y1 +
                 z1 * z1);
    double d2 = ((t[1][0] - px) * (t[1][0] - px) +
                 (t[1][1] - py) * (t[1][1] - py) +
                 (t[1][2] - pz) * (t[1][2] - pz));
    if(d2 < d1)
    {
        d1 = d2;
        X = t[1][0];
        Y = t[1][1];
        Z = t[1][2];
    }
    d2 = ((t[2][0] - px) * (t[2][0] - px) +
          (t[2][1] - py) * (t[2][1] - py) +
          (t[2][2] - pz) * (t[2][2] - pz));
    if(d2 < d1)
    {
        d1 = d2;
        X = t[2][0];
        Y = t[2][1];
        Z = t[2][2];
    }


    const double d = x3 * x3 * (y2 * y2 + z2 * z2) + (y3 * z2 - y2 * z3) * (y3 * z2 - y2 * z3) -
                     2 * x2 * x3 * (y2 * y3 + z2 * z3) + x2 * x2 * (y3 * y3 + z3 * z3);

    a = -((x3 * x3 * (y1 * y2 + z1 * z2) + (y3 * z1 - y1 * z3) * (y3 * z2 - y2 * z3) -
           x3 * (x2 * y1 * y3 + x1 * y2 * y3 + x2 * z1 * z3 + x1 * z2 * z3) +
           x1 * x2 * (y3 * y3 + z3 * z3)) / d);

    if(a >= 0 && a <= 1)
    {
        b = -((x1 * x3 * (y2 * y2 + z2 * z2) + (y2 * z1 - y1 * z2) * (-y3 * z2 + y2 * z3) +
               x2 * x2 * (y1 * y3 + z1 * z3) - x2 * (x3 * y1 * y2 + x1 * y2 * y3 + x3 * z1 * z2 +
                                                     x1 * z2 * z3)) / d);
        if(b >= 0 && (a + b) <= 1)
        {

            X = (x1 + x2 * a + x3 * b);
            Y = (y1 + y2 * a + y3 * b);
            Z = (z1 + z2 * a + z3 * b);
            return (X > -0.5) && (X < 0.5) && (Y > -0.5) && (Y < 0.5) && (Z > -0.5) && (Z < 0.5);
        }
    }

    a = -((x1 * x2 + y1 * y2 + z1 * z2) / (x2 * x2 + y2 * y2 + z2 * z2));
    if(a >= 0 && a <= 1)
    {
        d2 =    (x1 + x2 * a) * (x1 + x2 * a) +
                (y1 + y2 * a) * (y1 + y2 * a) +
                (z1 + z2 * a) * (z1 + z2 * a);
        if(d2 < d1)
        {
            b = 0;
            X = (x1 + x2 * a + x3 * b);
            Y = (y1 + y2 * a + y3 * b);
            Z = (z1 + z2 * a + z3 * b);
            d1 = d2;
        }
    }

    b = -((x1 * x3 + y1 * y3 + z1 * z3) / (x3 * x3 + y3 * y3 + z3 * z3));
    if(b >= 0 && b <= 1)
    {
        d2 =    (x1 + x3 * b) * (x1 + x3 * b) +
                (y1 + y3 * b) * (y1 + y3 * b) +
                (z1 + z3 * b) * (z1 + z3 * b);
        if(d2 < d1)
        {
            a = 0;
            X = (x1 + x2 * a + x3 * b);
            Y = (y1 + y2 * a + y3 * b);
            Z = (z1 + z2 * a + z3 * b);
            d1 = d2;
        }
    }

    a = -(((x2 - x3) * (x1 + x3) + (y2 - y3) * (y1 + y3) + (z2 - z3) *
                                                           (z1 + z3)) / ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3)));
    if(a >= 0 && a <= 1)
    {
        b = 1 - a;
        d2 =    (x1 + x2 * a + x3 * b) * (x1 + x2 * a + x3 * b) +
                (y1 + y2 * a + y3 * b) * (y1 + y2 * a + y3 * b) +
                (z1 + z2 * a + z3 * b) * (z1 + z2 * a + z3 * b);
        if(d2 < d1)
        {
            X = (x1 + x2 * a + x3 * b);
            Y = (y1 + y2 * a + y3 * b);
            Z = (z1 + z2 * a + z3 * b);
        }
    }
    return (X > -0.5) && (X < 0.5) && (Y > -0.5) && (Y < 0.5) && (Z > -0.5) && (Z < 0.5);
}

std::vector<std::array<size_t, 3>> find_seeds_fim(triangle *mesh, unsigned long X, unsigned long Y, unsigned long Z, size_t mesh_size) {
    auto grid_points = new bool[X * Y * Z];
    for(int i = 0; i < mesh_size; ++i)
    {
        find_triangle_grid_points(grid_points, mesh[i], X, Y, Z);
    }
    std::vector<std::array<size_t, 3>> result;
    for(unsigned long long z = 0; z < Z; ++z)
    {
        for(unsigned long long y = 0; y < Y; ++y)
        {
            for(unsigned long long x = 0; x < X; ++x)
            {
                if(grid_points[z * (X * Y) + y * X + x])
                {
                    result.emplace_back(std::array<size_t, 3>({x, y, z}));
                }
            }
        }
    }
    delete[] grid_points;
    return result;
}

bool outside_mesh(triangle *mesh, int px, int py, int pz, size_t mesh_size, int Y) {
    size_t intersected = 0;
    for(size_t i = 0; i < mesh_size; ++i)
    {
        const double    x1 = mesh[i][0][0] - px, y1 = mesh[i][0][1] - py, z1 = mesh[i][0][2] - pz,
                x2 = mesh[i][1][0] - mesh[i][0][0], y2 = mesh[i][1][1] - mesh[i][0][1], z2 = mesh[i][1][2] - mesh[i][0][2],
                x3 = mesh[i][2][0] - mesh[i][0][0], y3 = mesh[i][2][1] - mesh[i][0][1], z3 = mesh[i][2][2] - mesh[i][0][2];

        const double    k = 1 - (2.0 * py / Y);
        const double    a_inter2    = (k * x3 * z1 - y3 * z1 - k * x1 * z3 + y1 * z3) /
                                        (y3 * z2 - k * x3 * z2 + k * x2 * z3 - y2 * z3),
                        b_inter2    = -1 * (k * x2 * z1 - y2 * z1 - k * x1 * z2 + y1 * z2) /
                                (y3 * z2 - k * x3 * z2 + k * x2 * z3 - y2 * z3),
                        x_inter     = (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3)/
                                        (y3 * z2 - k * x3 * z2 + k * x2 * z3 - y2 * z3);
        if(a_inter2 > 0 && b_inter2 > 0 && (a_inter2 + b_inter2) < 1 && x_inter > 0) ++intersected;
    }
    return intersected % 2 == 0;
}
