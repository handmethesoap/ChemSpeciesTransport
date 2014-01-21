#include <Array.hh>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <utility>

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

Array::Array( int xSize, int ySize, int zSize)
{
    size_ = zSize*ySize*xSize;
    domain_ = new T[size_];
    xSize_ = xSize;
    ySize_ = ySize;
    zSize_ = zSize;
}

Array::Array(const Array &arr) : xSize_(arr.getSize(0)), ySize_(arr.getSize(1)), zSize_(arr.getSize(2)), size_(arr.getSize()), domain_(new T[arr.getSize()]) {
    std::copy(arr.domain_, arr.domain_ + arr.size_, domain_);
}

Array& Array::operator= (const Array& arr) {
    Array cpy(arr);
    std::swap(xSize_, cpy.xSize_);
    std::swap(ySize_, cpy.ySize_);
    std::swap(zSize_, cpy.zSize_);
    std::swap(size_, cpy.size_);
    std::swap(domain_, cpy.domain_);
    return *this;
}

Array::~Array() {
    delete[] domain_;
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( T value )
{
    std::fill(domain_,domain_ + size_, value);
}

std::pair<T,T> Array::minmax() {
    T minimum = 1/0.0;
    T maximum = -1/0.0;
    for (int k = 0; k < zSize_; k++) {
        for (int j = 0; j < ySize_; j++) {
            for (int i = 0; i < xSize_; i++) {
                maximum = std::max((*this)(i,j,k), maximum);
                minimum = std::min((*this)(i,j,k), minimum);
            }
        }
    }
    return std::make_pair(minimum,maximum);    
}

//initialize the whole array with a function
void Array::fill (std::function<T(int,int,int)> f) {
    for (int k = 0; k < zSize_; k++) {
        for (int j = 0; j < ySize_; j++) {
            for (int i = 0; i < xSize_; i++) {
                (*this)(i,j,k) = f(i,j,k);
            }
        }
    }
}

// Print the whole array (for debugging purposes)
void Array::print()
{
    // For 2D Arrays the positive x-coordinate goes to the right
    //                   positive y-coordinate goes upwards
    //      -> the line with highest y-value should be printed first
    for (int k = getSize(2)-1; k >= 0; k--) {
        std::cout << "z: " << k << std::endl;
        for (int j = getSize(1)-1; j >= 0; j--) {
            for (int i = 0; i < getSize(0); i++) {
                std::cout << std::fixed << (*this)(i,j,k) << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

int Array::getSize( int dimension ) const
{
    return (dimension == 0) ? xSize_ : ((dimension == 1) ? ySize_ : ((dimension == 2) ? zSize_ : -1));
}

//return total size of the array
int Array::getSize() const
{
    return size_;
}

bool Array::checkNaN() const {
   for (int k = 0; k < zSize_; k++) {
        for (int j = 0; j < ySize_; j++) {
            for (int i = 0; i < xSize_; i++) {
                T temp = (*this)(i,j,k);
                if (temp != temp) return false;
            }
        }
    }
    return true;
}
