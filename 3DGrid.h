#pragma once
#include <memory>
#include <stdexcept>
    

struct PointXYZ {
    float x, y, z;
    PointXYZ(float _x = 0.0, float _y = 0.0, float _z = 0.0) : x(_x), y(_y), z(_z) {}
};

class Array3D {
public:
    Array3D(int d, int h, int w) : depth(d), height(h), width(w)
    {
        data = std::make_unique<PointXYZ[]>(d * h * w);
    }

    void setPoint(int d, int h, int w, const PointXYZ& point) 
    {
        if (d >= 0 && d < depth && h >= 0 && h < height && w >= 0 && w < width) 
        {
            data[d * height * width + h * width + w] = point;
        }
        else
        {
            throw std::runtime_error("Index out of bounds!");
        }
    }

    PointXYZ getPoint(int d, int h, int w) const {
        if (d >= 0 && d < depth && h >= 0 && h < height && w >= 0 && w < width) 
        {
            return data[d * height * width + h * width + w];
        }
        throw std::runtime_error("Index out of bounds!");
    }

private:
    std::unique_ptr<PointXYZ[]> data;
    int depth, height, width;

};