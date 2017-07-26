// Mandelbrot.h
// Created by Laura Watkins.

#include <math.h>
#include "Grid.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>

#ifndef MANDELBROT
#define MANDELBROT

class Mandelbrot
{

public:

    Mandelbrot(double rad_max, int iter_max) :
        rad_max_(rad_max), iter_max_(iter_max) {}

    KOKKOS_INLINE_FUNCTION
    double calculate_mandelbrot(double x0, double y0)const {
        int iter = 0;
        double rad_z2 = 0.0;
        double x = 0.0;
        double y = 0.0;
        while ( rad_z2 < pow(rad_max_, 2) && iter < iter_max_ ) {
            double xtemp = x * x - y * y + x0 ;
            y = 2*x*y + y0 ;
            x = xtemp ;
            rad_z2 = x*x + y*y ;
            iter = iter + 1 ;
        }
        if ( iter < iter_max_) {
            return 0.;
        }
        else {
            return 255.;
        }
    }

    //Grayscale mandelbrot
    KOKKOS_INLINE_FUNCTION
    double calculate_mandelbrot_gs(double x0, double y0, Grid grid)const {
        int iter = 0;
        double rad_z2 = 0.0;
        double x = 0.0;
        double y = 0.0;
        double dx = 0.0;
        double dy = 0.0;
        while ( rad_z2 < pow(rad_max_,2) && iter < iter_max_ ) {
            dx = 2*x*dx + 1 ;
            dy = 2*y*dy ;
            double xtemp = x*x - y*y + x0 ;
            y = 2*x*y + y0 ;
            x = xtemp ;
            rad_z2 = x*x + y*y ;
            iter = iter + 1 ;
        }
        double rad_z = sqrt(rad_z2) ;
        double dz = sqrt(dx*dx + dy*dy) ;
        double distance = 2*log(rad_z)*rad_z ;
        double gs_value = color(distance, grid) ;
        return gs_value ;
    }

    // Get grayscale color between 0 and 255
    double color(double distance, Grid grid)const
    {
        double half_pixel_size = 0.5 * grid.pixel_size;
    
        if( distance < half_pixel_size )
        {
            return pow(distance / half_pixel_size, 1. / 3.) * 255.;
        }
        else
            return 255.;
    }
    
private:

    //KOKKOS_INLINE_FUNCTION
    //double color(double distance, Grid grid);

    double rad_max_;

    int iter_max_;

};

#endif
