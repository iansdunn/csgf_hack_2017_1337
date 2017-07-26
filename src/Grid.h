// Grid.h
// Created by Ian Dunn.

#include <utility>
#include <math.h>

class Grid
{

public:

    Grid(double center_x, double center_y,
         double length_x, double length_y, 
         int pixel_count_x) :
        center_x_(center_x), center_y_(center_y),
        length_x_(length_x), length_y_(length_y),
        pixel_count_x_(pixel_count_x)
    {
        setup();
    }

    std::pair<double, double> calculate_xy(int i);

    int num_pixels;

private:

    void setup();

    int pixel_count_x_;
    int pixel_count_y_;
    
    double center_x_;
    double center_y_;

    double length_x_;
    double length_y_;

    double pixel_size_;

    double min_x_;
    double max_y_;

};
