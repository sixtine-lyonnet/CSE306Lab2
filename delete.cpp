#include <flann/flann.hpp>
#include <iostream>

int main() {
    float dataset[4][2] = {{1,2},{3,4},{5,6},{7,8}};
    flann::Matrix<float> data(&dataset[0][0], 4, 2);
    flann::Index<flann::L2<float>> index(data, flann::KDTreeIndexParams(1));
    index.buildIndex();
    std::cout << "FLANN index built." << std::endl;
    return 0;
}