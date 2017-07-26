// Mandelbrot.h
// Created by Laura Watkins.

#include <math.h>
#include "Grid.h"

class Mandelbrot
{

public:

    Mandelbrot(double rad_max, int iter_max) :
        rad_max_(rad_max), iter_max_(iter_max) {}
    int calculate_mandelbrot(double x0, double y0);

private:

    double color(double distance, Grid* grid);

    double rad_max_;

    int iter_max_;

};
