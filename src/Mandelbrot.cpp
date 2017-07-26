// Mandelbrot.cpp
// Created by Laura Watkins.

include "Mandelbrot.h"

int Mandelbrot::calculate_mandelbrot(double x0, double y0) {
    int iter = 0;
    double rad_z = 0.0;
    double x = 0.0;
    double y = 0.0;
    while ( rad_z < rad_max**2 && iter < iter_max ) {
        xtemp = x*x - y*y + x0 ;
        y = 2*x*y + y0 ;
        x = xtemp ;
        rad_z = x*x + y*y ;
        iter = iter + 1 ;
    }
    if ( iter < iter_max) {
        return 0;
    }
    else {
        return 1;
    }
}
