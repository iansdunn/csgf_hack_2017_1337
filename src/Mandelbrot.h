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
};

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

    KOKKOS_INLINE_FUNCTION
    mytriple hsv_to_rgb(mytriple hsv)const
    {
        int h = hsv.a;
        int s = hsv.b;
        int v = hsv.c;

        int c = v * s;
        int hp = floor(h / 60.);

        int x = c * (hp % 2);

        mytriple rgb;

        switch (hp)
        {
            case 0:
                rgb = {c, x, 0};
            case 1:
                rgb = {x, c, 0};
            case 2:
                rgb = {0, c, x};
            case 3:
                rgb = {0, x, c};
            case 4:
                rgb = {x, 0, c};
            case 5:
                rgb = {c, 0, x};
            default:        
                rgb = {0, 0, 0};
}

        int m = v - c;

        rgb.a += m;
        rgb.b += m;
        rgb.c += m;

        return rgb;
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
        double gs_value = graysc(distance, grid) ;
        return gs_value ;
    }

    //COLOR mandelbrot
    KOKKOS_INLINE_FUNCTION
    mytriple calculate_mandelbrot_color(double x0, double y0, Grid grid)const {
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
        return color(distance, iter, grid) ;
    }

    // Get grayscale color between 0 and 255
    KOKKOS_INLINE_FUNCTION
    double graysc(double distance, Grid grid)const
    {
        double half_pixel_size = 0.5 * grid.pixel_size;
    
        if( distance < half_pixel_size )
        {
            return pow(distance / half_pixel_size, 1. / 3.) * 255.;
        }
        else
            return 255.;
    }

    KOKKOS_INLINE_FUNCTION
    mytriple color(double distance, int iterations, Grid grid)const
    {
        double half_pixel_size = 0.5 * grid.pixel_size;

        double value = 0;

        if( iterations > iter_max_) {
            mytriple hsv = {0, 0, 0} ;
            return hsv_to_rgb(hsv) ;
        } else if( distance < half_pixel_size ) {
            value = pow(distance/ half_pixel_size, 1./3.) ;
        } else {
            value = 1.0 ;
        }
       
        double saturation = 0.7 ; 
        double hue = log(iterations)/log(iter_max_) ;
        hue = hue - floor(hue) ;

        mytriple hsv = {hue, saturation, value} ;
        return hsv_to_rgb(hsv) ;
    }
    
private:

    double rad_max_;

    int iter_max_;

};

#endif
