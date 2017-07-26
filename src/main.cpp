#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>
#include "mpi.h"
#include <assert.h>
#include <limits>
#include "Grid.h" 
#include "Mandelbrot.h"
#include <fstream>
#include <bitset>

using namespace std;

class MyFunctor { 
public:
  MyFunctor(Kokkos::View<mytriple*> a, Grid grid, Mandelbrot mb) : a_(a), grid_(grid), mb_(mb) {}
  
  KOKKOS_INLINE_FUNCTION 
  void operator() (int i) const {  
    mypair xy = grid_.calculate_xy(i);                      
    a_(i) = mb_.calculate_mandelbrot_color(xy.x, xy.y, grid_);
  }

private:
  Kokkos::View<mytriple*> a_;

  Grid grid_;

  Mandelbrot mb_;
};

int main(int argc, char **argv) {
  // Initialize MPI before Kokkos
  //MPI_Init(&argc, &argv);

  // Initialize Kokkos
  Kokkos::initialize(argc, argv);
  
  // Set up grid
  const double center_x = -0.75;
  const double center_y = 0.00;
  const double length_x = 2.75;
  const double length_y = 2.0;
  const int pixel_count_x = 1000;
  Grid grid(center_x, center_y, length_x, length_y, pixel_count_x);
  
  // Set up Mandelbrot object
  const double rad_max = 2.0;
  const int iter_max = 1000;
  Mandelbrot mandelbrot(rad_max, iter_max);

  // Allocate our arrays
  Kokkos::View<mytriple*> a("a", grid.num_pixels);

  // Create host mirror of a
  auto  a_mirror = Kokkos::create_mirror_view(a);

  // parallel loop
  MyFunctor body = MyFunctor(a, grid, mandelbrot);
  Kokkos::parallel_for(grid.num_pixels, body); 

  // Update the mirror
  deep_copy(a_mirror, a);
  
  FILE* fid = fopen("mandelbrot.ppm", "wb");
  fprintf(fid, "P5\n");
  fprintf(fid, "%i %i\n", grid.pixel_count_x, grid.pixel_count_y);
  fprintf(fid, "%i\n", 255);
  for (auto i = 0; i < grid.num_pixels; i++) {
    unsigned char tmpa = a_mirror(i).a;
    unsigned char tmpb = a_mirror(i).b;
    unsigned char tmpc = a_mirror(i).c;
    fwrite(&tmpa, sizeof(unsigned char), 1, fid);
    fwrite(&tmpb, sizeof(unsigned char), 1, fid);
    fwrite(&tmpc, sizeof(unsigned char), 1, fid);
  }
  fclose(fid);
  
  cout << "Done!" << endl;

  Kokkos::finalize();

  //MPI_Finalize();

  return 0;
}
