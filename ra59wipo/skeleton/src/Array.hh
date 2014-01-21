#ifndef ARRAY_HH
#define ARRAY_HH


#include <Types.hh>
#include <Debug.hh>
#include <functional>
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <utility>
//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
template <class T>
class Array
{
public:
    // Constructors for 1D,2D and 3D
    Array( int xSize = 1, int ySize = 1, int zSize = 1 );

    // Depending on your implementation you might need the following:
    ~Array();
    Array(const Array& s);
    Array& operator= (const Array& s);

    // Access Operators for 1D, 2D and 3D
    inline T & operator () ( int i, int j = 0, int k = 0 );

    // for const Arrays the following access operators are required
    inline const T & operator () ( int i, int j = 0, int k = 0 ) const;

    // initialize the whole array with a constant value
    void fill( T value );

    // initialize the whole array with a lambda function
    void fill (std::function<T(int,int,int)> f);

    // return total size of the array
    int getSize() const;

    // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
    // other dimension values are not allowed
    int getSize(int dimension ) const;

    // Print the whole array ( for debugging purposes )
    void print();
    
    bool checkNaN() const; 
    
    std::pair<T,T> minmax();

private:
    int xSize_;
    int ySize_;
    int zSize_;
    int size_;
    T *domain_;
};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================

// Operator() 1D,2D,3D
template <class T>
inline T& Array<T>::operator ()(int i, int j, int k)
{
    ASSERT_MSG(((i >= 0) && (i < xSize_) && (j >= 0) && (j < ySize_) && (k >= 0) && (k < zSize_)), "Index out of bounds.")
    return domain_[i + j * xSize_ + k * xSize_ * ySize_];
}

template <class T>
inline const T& Array<T>::operator ()(int i, int j, int k) const
{
    ASSERT_MSG(((i >= 0) && (i < xSize_) && (j >= 0) && (j < ySize_) && (k >= 0) && (k < zSize_)), "Index out of bounds.")
    return domain_[i + j * xSize_ + k * xSize_ * ySize_];
}

template <class T>
Array<T>::Array( int xSize, int ySize, int zSize)
{
    size_ = zSize*ySize*xSize;
    domain_ = new T[size_];
    xSize_ = xSize;
    ySize_ = ySize;
    zSize_ = zSize;
}

template <class T>
Array<T>::Array(const Array &arr) : xSize_(arr.getSize(0)), ySize_(arr.getSize(1)), zSize_(arr.getSize(2)), size_(arr.getSize()), domain_(new T[arr.getSize()]) {
    std::copy(arr.domain_, arr.domain_ + arr.size_, domain_);
}

template <class T>
Array<T>& Array<T>::operator= (const Array<T>& arr) {
    Array cpy(arr);
    std::swap(xSize_, cpy.xSize_);
    std::swap(ySize_, cpy.ySize_);
    std::swap(zSize_, cpy.zSize_);
    std::swap(size_, cpy.size_);
    std::swap(domain_, cpy.domain_);
    return *this;
}

template <class T>
Array<T>::~Array() {
    delete[] domain_;
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value


template <class T>
void Array<T>::fill( T value )
{
    std::fill(domain_,domain_ + size_, value);
}

template <class T>
std::pair<T,T> Array<T>::minmax() {
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
template <class T>
void Array<T>::fill (std::function<T(int,int,int)> f) {
    for (int k = 0; k < zSize_; k++) {
        for (int j = 0; j < ySize_; j++) {
            for (int i = 0; i < xSize_; i++) {
                (*this)(i,j,k) = f(i,j,k);
            }
        }
    }
}

// Print the whole array (for debugging purposes)
template <class T>
void Array<T>::print()
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

template <class T>
int Array<T>::getSize( int dimension ) const
{
    return (dimension == 0) ? xSize_ : ((dimension == 1) ? ySize_ : ((dimension == 2) ? zSize_ : -1));
}

//return total size of the array
template <class T>
int Array<T>::getSize() const
{
    return size_;
}

template <class T>
bool Array<T>::checkNaN() const {
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

#endif //ARRAY_HH

