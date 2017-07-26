// Grid.h
// Created by Ian Dunn.

class Grid
{

public:

    Grid();
    ~Grid();

    void setup();

private:

    double center_x_;
    double center_y_;

    double length_x_;
    double length_y_;

    int pixel_count_x_;
    int pixel_count_y_;
    
    double pixel_size_;

    double min_x_;
    double max_y_;

}
