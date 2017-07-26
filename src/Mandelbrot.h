// Mandelbrot.h
// Created by Laura Watkins.

#include <math.h>
#include "Grid.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>

#ifndef MANDELBROT
#define MANDELBROT

struct mytriple
{
    int a;
    int b;
    int c;
}

class Mandelbrot
{

public:

    Mandelbrot(double rad_max, int iter_max) :
        rad_max_(rad_max), iter_max_(iter_max) {}

    KOKKOS_INLINE_FUNCTION
    int calculate_mandelbrot(double x0, double y0)const {
        int iter = 0;
        double rad_z = 0.0;
        double x = 0.0;
        double y = 0.0;
        while ( rad_z < pow(rad_max_, 2) && iter < iter_max_ ) {
            double xtemp = x * x - y * y + x0 ;
            y = 2*x*y + y0 ;
            x = xtemp ;
            rad_z = x*x + y*y ;
            iter = iter + 1 ;
        }
        if ( iter < iter_max_) {
            return 0;
        }
        else {
            return 255;
        }
    }

    KOKKOS_INLINE_FUNCTION
    mytriple hsv_to_rgb(mytriple hsv)
    {
        mytriple rgb;
        rgb.a = 0;
        rgb.b = 0;
        rgb.c = 0;

        return rgb.c;
    }
    
private:

    //KOKKOS_INLINE_FUNCTION
    //double color(double distance, Grid grid);

    double rad_max_;

    int iter_max_;

};

#endif
