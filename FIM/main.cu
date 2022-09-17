#include <StructuredEikonal.h>
#include <happly.h>
#include "../mesh_functions.h"
#include "../FMM/mesh_functions_fmm.h"
#include <chrono>
#include <cmath>

int main()
{
    happly::PLYData plyIn(R"(C:\Users\Maks\CLionProjects\Agate3DEvolver\AgateContour_0209.ply)");
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
    /*const int min_x = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[0] < b[0]);}))[0]);
    const int min_y = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[1] < b[1]);}))[1]);
    const int min_z = std::floor((*std::min_element(vPos.begin(), vPos.end(),
                                                    [](const std::array<double, 3>& a, const std::array<double, 3>& b)
                                                    {return (a[2] < b[2]);}))[2]);
                                                    */
    const int min_x = 0, min_y = 0, min_z = 0;
    const size_t X = max_x - min_x, Y = max_y - min_y, Z = max_z - min_z;
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

    auto boundary_indices = find_seeds_fim(mesh, X, Y, Z, fInd.size());

    std::cout << fInd.size() << " triangles loaded from mesh! Starting calculations... \n";

    StructuredEikonal data(false);
    data.setDims(360, 360, 360);
    data.setSeeds(boundary_indices);
    data.solveEikonal();
    auto arrival_times = data.answer_;

    for(int z = 0; z < Z; ++z)
    {
        for(int y = 0; y < Y; ++y)
        {
            for(int x = 0; x < X; ++x)
            {
                if(outside_mesh(mesh, x, y, z, fInd.size(), Y))
                {
                    arrival_times[x][y][z] = -std::fabs(arrival_times[x][y][z]);
                }
                else
                {
                    arrival_times[x][y][z] = std::fabs(arrival_times[x][y][z]);
                }
            }
        }
    }

    delete[] mesh;

    /*
    void StructuredEikonal::setDims(size_t w, size_t h, size_t d);  //set the volume dimensions
    void StructuredEikonal::setMapType(size_t t); //pre-generated speed functions (sphere or egg-carton)
    void StructuredEikonal::setItersPerBlock(size_t t); //set the iterations per block
    void StructuredEikonal::setSpeeds(std::vector<std::vector<std::vector<double> > > speed); //set the voxel speeds
    void StructuredEikonal::setSeeds(std::vector<std::array<size_t, 3> > seeds); //set list of seed voxels
     */
    auto stop = std::chrono::high_resolution_clock::now();


    std::ofstream logfile("log.txt", std::ios::app);
    logfile << "FIM " << std::chrono::duration<double, std::ratio<3600>>(stop - start).count() << "\n";
    logfile.close();

    start = std::chrono::high_resolution_clock::now();
    std::ofstream out("AgateDistanceMap_fim3.txt", std::ios::out);
    for(int z = 0; z < 360; ++z)
    {
        for(int y = 0; y < 360; ++y)
        {
            for(int x = 0; x < 360; ++x)
            {
                out << x << ' ' << y << ' ' << z << ' ' << arrival_times[x][y][z] << '\n';
            }
        }
    }
    out.close();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Writing to file took: " << std::chrono::duration<double>(stop - start).count() <<
              " seconds\n";
}