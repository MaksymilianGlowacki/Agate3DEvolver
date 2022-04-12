#include "mesh_functions_gpu.cuh"

__device__
double distance_to_triangle(triangle t, double px, double py, double pz, size_t& inter)
{
    double a, b;
    double d1 = ((t[0][0] - px) * (t[0][0] - px) +
                 (t[0][1] - py) * (t[0][1] - py) +
                 (t[0][2] - pz) * (t[0][2] - pz));
    double d2 = ((t[1][0] - px) * (t[1][0] - px) +
                 (t[1][1] - py) * (t[1][1] - py) +
                 (t[1][2] - pz) * (t[1][2] - pz));
    if(d2 < d1) d1 = d2;
    d2 = ((t[2][0] - px) * (t[2][0] - px) +
          (t[2][1] - py) * (t[2][1] - py) +
          (t[2][2] - pz) * (t[2][2] - pz));
    if(d2 < d1) d1 = d2;

    const double    x1 = t[0][0] - px, y1 = t[0][1] - py, z1 = t[0][2] - pz,
                    x2 = t[1][0] - t[0][0], y2 = t[1][1] - t[0][1], z2 = t[1][2] - t[0][2],
                    x3 = t[2][0] - t[0][0], y3 = t[2][1] - t[0][1], z3 = t[2][2] - t[0][2];

    const double    a_inter = -(x3 * y1 - x1 * y3) / (x3 * y2 - x2 * y3),
                    b_inter = (x2 * y1 - x1 * y2) / (x3 * y2 - x2 * y3),
                    x_inter = (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3)/
                                (x3 * y2 - x2 * y3);
     if(a_inter > 0 && b_inter > 0 && (a_inter + b_inter) < 1 && x_inter > 0) ++inter;

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
            return std::sqrt((x1 + x2 * a + x3 * b) * (x1 + x2 * a + x3 * b) +
                             (y1 + y2 * a + y3 * b) * (y1 + y2 * a + y3 * b) +
                             (z1 + z2 * a + z3 * b) * (z1 + z2 * a + z3 * b));
        }
    }

    a = -((x1 * x2 + y1 * y2 + z1 * z2) / (x2 * x2 + y2 * y2 + z2 * z2));
    if(a >= 0 && a <= 1)
    {
        d2 =    (x1 + x2 * a) * (x1 + x2 * a) +
                (y1 + y2 * a) * (y1 + y2 * a) +
                (z1 + z2 * a) * (z1 + z2 * a);
        if(d2 < d1) d1 = d2;
    }

    b = -((x1 * x3 + y1 * y3 + z1 * z3) / (x3 * x3 + y3 * y3 + z3 * z3));
    if(b >= 0 && b <= 1)
    {
        d2 =    (x1 + x3 * b) * (x1 + x3 * b) +
                (y1 + y3 * b) * (y1 + y3 * b) +
                (z1 + z3 * b) * (z1 + z3 * b);
        if(d2 < d1) d1 = d2;
    }

    a = -(((x2 - x3) * (x1 + x3) + (y2 - y3) * (y1 + y3) + (z2 - z3) *
            (z1 + z3)) / ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3)));
    if(a >= 0 && a <= 1)
    {
        b = 1 - a;
        d2 =    (x1 + x2 * a + x3 * b) * (x1 + x2 * a + x3 * b) +
                (y1 + y2 * a + y3 * b) * (y1 + y2 * a + y3 * b) +
                (z1 + z2 * a + z3 * b) * (z1 + z2 * a + z3 * b);
        if(d2 < d1) d1 = d2;
    }
    return std::sqrt(d1);
}

__device__
double distance_to_mesh(triangle *mesh, unsigned int x, unsigned int y, unsigned int z, size_t mesh_size)
{
    size_t intersected = 0;
    double d_min = distance_to_triangle(mesh[0], x, y, z, intersected), d;
    for(size_t i = 1; i < mesh_size; ++i)
    {
        d = distance_to_triangle(mesh[i], x, y, z, intersected);
        if(d < d_min) d_min = d;
    }
    return (intersected % 2 == 0) ? d_min : -d_min;
}