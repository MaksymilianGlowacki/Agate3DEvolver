#include "mesh_functions_gpu.cuh"
#include <happly.h> //edited lines 954 & 997: dynamic_cast converted to static_cast (avoid compilation error)
#include <algorithm>
#include <chrono>

__global__
void calculate_distances(const int X, const int Y, size_t num_tri, triangle *mesh, double *results)
{
    unsigned int z = threadIdx.x;
    unsigned int y = blockIdx.y;
    unsigned int x = blockIdx.x;
    results[z * (X * Y) + y * X + x] = distance_to_mesh(mesh, x, y, z, num_tri);
}

int main()
{
    happly::PLYData plyIn(R"(C:\Users\Maks\CLionProjects\Agate3DEvolver\AgateContour_1008.ply)");
    std::vector<std::array<double, 3>> vPos = plyIn.getVertexPositions();
    std::vector<std::vector<size_t>> fInd = plyIn.getFaceIndices<size_t>();

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
    const int X = max_x - min_x, Y = max_y - min_y, Z = max_z - min_z;
    triangle *mesh;
    double *results;
    cudaMallocManaged(&mesh, fInd.size()*sizeof(triangle));
    cudaMallocManaged(&results, (X * Y * Z)*sizeof(double));

    for(size_t i = 0; i < fInd.size(); ++i)
    {
        mesh[i][0][0] = vPos[fInd[i][0]][0] - min_x;
        mesh[i][0][1] = vPos[fInd[i][0]][1] - min_y;
        mesh[i][0][2] = vPos[fInd[i][0]][2] - min_z;

        mesh[i][1][0] = vPos[fInd[i][1]][0] - min_x;
        mesh[i][1][1] = vPos[fInd[i][1]][1] - min_y;
        mesh[i][1][2] = vPos[fInd[i][1]][2] - min_z;

        mesh[i][2][0] = vPos[fInd[i][2]][0] - min_x;
        mesh[i][2][1] = vPos[fInd[i][2]][1] - min_y;
        mesh[i][2][2] = vPos[fInd[i][2]][2] - min_z;
    }
    std::cout << fInd.size() << " triangles loaded from mesh! Starting calculations... \n";
    auto start = std::chrono::high_resolution_clock::now();

    calculate_distances<<<dim3(X, Y), Z>>>(X, Y, fInd.size(), mesh, results);

    cudaDeviceSynchronize();
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Calculations took: " << std::chrono::duration<double, std::ratio<3600>>(stop - start).count() <<
              " hours\nSaving results to .txt file...\n";
    start = std::chrono::high_resolution_clock::now();
    std::ofstream out("AgateDistanceMap3.txt", std::ios::out);
    for(int z = 0; z < Z; ++z)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int x = 0; x < X; ++x)
            {
                out << x << ' ' << y << ' ' << z << ' ' << results[z * (X * Y) + y * X + x] << '\n';
            }
        }
    }
    cudaFree(mesh);
    cudaFree(results);
    out.close();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Writing to file took: " << std::chrono::duration<double>(stop - start).count() <<
              " seconds\n";
    return 0;
}
