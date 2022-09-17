#include <fast_marching_method.hpp>
#include <happly.h>
#include <thread_pool.hpp>
#include "../mesh_functions.h"
#include "mesh_functions_fmm.h"

namespace fmm = thinks::fast_marching_method;

int main(int argc, char* argv[])
{
    int num_threads = 16, factor = 20;
    if(argc > 2)
    {
        num_threads = std::atol(argv[1]);
        factor = std::atol(argv[2]);
    }

    happly::PLYData plyIn(R"(C:\Users\Maks\CLionProjects\Agate3DEvolver\AgateContour_0209.ply)");
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();

    /*
    const int max_x = std::ceil((*std::max_element(vPos.begin(), vPos.end(),
                                                   [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                   {return (a[0] < b[0]);}))[0]);
    const int max_y = std::ceil((*std::max_element(vPos.begin(), vPos.end(),
                                                   [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                   {return (a[1] < b[1]);}))[1]);
    const int max_z = std::ceil((*std::max_element(vPos.begin(), vPos.end(),
                                                   [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                   {return (a[2] < b[2]);}))[2]);
    const int min_x = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[0] < b[0]);}))[0]);
    const int min_y = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[1] < b[1]);}))[1]);
    const int min_z = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[2] < b[2]);}))[2]);
    const unsigned long X = max_x - min_x, Y = max_y - min_y, Z = max_z - min_z;
     */
    const unsigned long X = 300, Y = 210, Z = 361, min_x = 0, min_y = 0, min_z = 0;
    auto *mesh = new triangle[fInd.size()];
    for(size_t i = 0; i < fInd.size(); ++i)
    {
        mesh[i][0][0] = vPos[fInd[i][0]][0] - min_x;
        mesh[i][0][1] = vPos[fInd[i][0]][1] - min_y;
        mesh[i][0][2] = vPos[fInd[i][0]][2] - min_z + 0.4;

        mesh[i][1][0] = vPos[fInd[i][1]][0] - min_x;
        mesh[i][1][1] = vPos[fInd[i][1]][1] - min_y;
        mesh[i][1][2] = vPos[fInd[i][1]][2] - min_z + 0.4;

        mesh[i][2][0] = vPos[fInd[i][2]][0] - min_x;
        mesh[i][2][1] = vPos[fInd[i][2]][1] - min_y;
        mesh[i][2][2] = vPos[fInd[i][2]][2] - min_z + 0.4;
    }
    auto start = std::chrono::high_resolution_clock::now();

    auto boundary_indices2 = find_seeds(mesh, X, Y, Z, fInd.size());
    std::vector<std::array<int32_t, 3>> boundary_indices;
    for(int i = 0; i < boundary_indices2.size(); i += factor)
    {
        boundary_indices.emplace_back(boundary_indices2[i]);
    }
    auto *boundary_distances2 = new std::future<double>[boundary_indices.size()];
    std::vector<float> boundary_distances(boundary_indices.size());
    //boundary_distances.reserve(boundary_indices.size());
    /*
    for(const auto & point : boundary_indices)
    {
        boundary_distances.emplace_back((float)distance_to_mesh(mesh, point[0], point[1], point[2], fInd.size(), Y));
    }
     */
    thread_pool pool(num_threads);
    for(int i = 0; i < boundary_indices.size(); ++i)
    {
        boundary_distances2[i] = pool.submit(distance_to_mesh, mesh, boundary_indices[i][0],
                                                             boundary_indices[i][1],
                                                             boundary_indices[i][2], fInd.size(), Y);
    }
    pool.wait_for_tasks();

    for(int i = 0; i < boundary_indices.size(); ++i)
    {
        boundary_distances[i] = (float)(boundary_distances2[i].get());
    }
    delete[] boundary_distances2;
    std::cout << boundary_indices.size() << ' ' << boundary_distances.size() << '\n';

    /*
    std::ofstream seeds("seeds2.txt", std::ios::out);
    for(int i = 0; i < boundary_indices.size(); ++i)
    {
        seeds << boundary_indices[i][0] << ' ' << boundary_indices[i][1] << ' ' << boundary_indices[i][2] << ' '
            << boundary_distances[i] << '\n';
    }
    seeds.close();
     */

    std::cout << fInd.size() << " triangles loaded from mesh! Starting calculations... \n";

    auto grid_size = std::array<size_t, 3>{{X, Y, Z}};
    auto grid_spacing = std::array<float, 3>{{1.0f, 1.0f, 1.0f}};
    auto uniform_speed = 1.0f;

    auto *arrival_times = new std::vector<float>;
    *arrival_times = fmm::SignedArrivalTime(
            grid_size,
            boundary_indices,
            boundary_distances,
            fmm::UniformSpeedEikonalSolver<float, 3>(grid_spacing, uniform_speed)
            //fmm::HighAccuracyUniformSpeedEikonalSolver<float, 3>(1)
                    );


    for(int z = 0; z < Z; ++z)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int x = 0; x < X; ++x)
            {
                if(outside_mesh(mesh, x, y, z, fInd.size(), Y))
                {
                    (*arrival_times)[z * (X * Y) + y * X + x] = -std::fabs((*arrival_times)[z * (X * Y) + y * X + x]);
                }
                else
                {
                    (*arrival_times)[z * (X * Y) + y * X + x] = std::fabs((*arrival_times)[z * (X * Y) + y * X + x]);
                }
            }
        }
    }

    delete[] mesh;

    auto stop = std::chrono::high_resolution_clock::now();


    std::ofstream logfile("log.txt", std::ios::app);
    logfile << "FMM " << std::chrono::duration<double, std::ratio<3600>>(stop - start).count() <<
            " " << num_threads << "\n";
    logfile.close();

    start = std::chrono::high_resolution_clock::now();
    std::ofstream out("AgateDistanceMap_fmm3.txt", std::ios::out);
    for(int z = 0; z < Z; ++z)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int x = 0; x < X; ++x)
            {
                out << x << ' ' << y << ' ' << z << ' ' << (*arrival_times)[z * (X * Y) + y * X + x] << '\n';
            }
        }
    }
    delete arrival_times;
    out.close();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Writing to file took: " << std::chrono::duration<double>(stop - start).count() <<
              " seconds\n";
}