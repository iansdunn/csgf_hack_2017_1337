// Grid.cpp
// Created by Ian Dunn.

std::pair<double, double> Grid::calculate_xy(int i)
{
    int pixel_x = i % pixel_count_x_;
    int pixel_y = round((i - pixel_count_x) / pixel_count_x_);

    double x = min_x_ + pixel_x * pixel_size_;
    double y = max_y_ - pixel_y * pixel_size_;

    std::pair<double, double> xy(x, y);

    return xy;
}

void Grid::setup()
{
    min_x_ = center_x_ - length_x_ / 2.0;
    max_y_ = center_y_ + length_y_ / 2.0;

    pixel_size_ = length_x_ / pixel_count_x_;

    pixel_count_y_ = ceil(length_y_ / pixel_size_);

    num_pixels = pixel_count_x_ * pixel_count_y_;
}
