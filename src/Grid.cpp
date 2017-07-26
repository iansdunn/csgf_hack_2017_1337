// Grid.cpp
// Created by Ian Dunn.

#include "Grid.h"

void Grid::setup()
{
    min_x_ = center_x_ - length_x_ / 2.0;
    max_y_ = center_y_ + length_y_ / 2.0;

    pixel_size = length_x_ / pixel_count_x_;

    pixel_count_y_ = ceil(length_y_ / pixel_size);

    num_pixels = pixel_count_x_ * pixel_count_y_;
}
