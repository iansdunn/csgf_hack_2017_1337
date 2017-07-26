#include <iostream>
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>
#include "mpi.h"
#include <assert.h>
#include <limits>
#include "grind.h" 

class MyFunctor { 
public:
  MyFunctor(Kokkos::View<double*> a) : a{a} {}
  
  KOKKOS_INLINE_FUNCTION 
  void operator() (int i) const {  
    std::pair<double, double> xy = grid.calculate_xy(i);                      
    a(i) = Mandelbrot.calculate_mandelbrot(xy.first, xy.second);
  }
};

private:
  Kokkos::View<double*> a;
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
  const int pixel_count_x = 10;
  Grid grid = Grid(center_x, center_y, length_x, length_y, pixel_count_x);
  
  // Set up Mandelbrot object
  const double rad_max = 2.0;
  const int iter_max = 1000;
  Mandelbrot mandelbrot = Mandelbrot(rad_max, iter_max);

  // Allocate our arrays
  Kokkos::View<double*> a("a", grid.num_pixels);

  // Create host mirror of a
  auto  a_mirror = Kokkos::create_mirror_view(a);

  // parallel loop
  MyFunctor body = MyFunctor(a);
  Kokkos::parallel_for(grid.num_pixels, body); 

  // Update the mirror
  deep_copy(a_mirror, a);

  std::cout << "Done!" << endl;

  Kokkos::finalize();

  //MPI_Finalize();

  return 0;
}
