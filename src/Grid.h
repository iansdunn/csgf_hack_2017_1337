// Grid.h
// Created by Ian Dunn.

#include <utility>
#include <math.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_View.hpp>

#ifndef GRID
#define GRID

struct mypair
{
    double x;
    double y;
};

class Grid
{

public:

    Grid(double center_x, double center_y,
         double length_x, double length_y, 
         int pixel_count_x) :
        center_x_(center_x), center_y_(center_y),
        length_x_(length_x), length_y_(length_y),
        pixel_count_x(pixel_count_x)
    {
        setup();
    }

    KOKKOS_INLINE_FUNCTION
    mypair calculate_xy(int i)const
    {
        int pixel_x = i % pixel_count_x;
        int pixel_y = (i - pixel_count_x) / pixel_count_x;
    
        double x = min_x_ + pixel_x * pixel_size;
        double y = max_y_ - pixel_y * pixel_size;
    
        mypair xy;
        xy.x = x;
        xy.y = y;
    
        return xy;
    }
    
    int num_pixels;

    double pixel_size;

    int pixel_count_x;
    int pixel_count_y;
    
private:

    void setup();

    double center_x_;
    double center_y_;

    double length_x_;
    double length_y_;

    double min_x_;
    double max_y_;

};

#endif //GRID
