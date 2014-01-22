#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


#include "Types.hh"
#include "Array.hh"
#include <FileReader.hh>
#include <cmath>
#include "GrayScaleImage.hh"

#include <sstream>

//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
   // Constructors to manually create staggered grid
   StaggeredGrid ( int iimax, int jjmax, real ddx, real ddy );

   // Constructor to create a staggered grid from a parsed configuration file
   StaggeredGrid ( const FileReader & configuration );
   
   static void registerModule( FileReader & configuration );

   // Getters / Setters for member variables
   Array<real> & p()    { return p_;    }
   Array<real> & rhs()  { return rhs_;  }
   Array<real> & u()    { return u_;  }
   Array<real> & v()    { return v_;  }
   Array<real> & f()    { return f_;  }
   Array<real> & g()    { return g_;  }
   Array<bool> & obstacleflags() { return obstacleflags_; }

   Array<real> & c(int i) { return c_[i]; }

   const Array<real> & p()   const { return p_;   }
   const Array<real> & rhs() const { return rhs_; }
   const Array<real> & u()   const { return u_;   }
   const Array<real> & v()   const { return v_;   }
   const Array<real> & f()   const { return f_;   }
   const Array<real> & g()   const { return g_;   }
   const Array<bool> & obstacleflags() const { return obstacleflags_; }

   const Array<real> & c(int i) const { return c_[i]; }
   
   // wrapped access
   inline real p(const int x, const int y, Direction dir);
   inline real u(const int x, const int y, Direction dir);
   inline real v(const int x, const int y, Direction dir);
   
   real dx() const { return dx_; }
   real dy() const { return dy_; }

   int imax() const { return imax_; }
   int jmax() const { return jmax_; }
  
    int numSpecies() const { return numSpecies_; }
    real lambda(int i) const { return lambda_[i]; }

   int xSize() const { return imax_; }
   int ySize() const { return jmax_; }
   
   inline bool isFluid (const int x, const int y);
   inline int getNumFluid();
   inline void setCellToObstacle(const int x, const int y);

   void createRectangle(int x1, int y1, int x2, int y2);
   void createCircle(int x, int y, int r);
   
   void saveObstacleFieldToPNG(std::string filename);
   void loadObstacleFieldFromPNG(std::string filename);
   
protected:
   Array<real> p_;   //< pressure field
   Array<real> rhs_; //< right hand side of the pressure equation
   Array<real> u_;   //< horizontal velocity component
   Array<real> v_;   //< vertical velocity compoent
   Array<real> f_;   //< f (helper array)
   Array<real> g_;   //< g (helper array)
   Array<bool> obstacleflags_;

   std::vector<real> lambda_;
   std::vector<Array<real> > c_;

   real dx_;   //< distance between two grid points in x direction
   real dy_;   //< distance between two grid points in y direction

   int imax_;
   int jmax_;
   
    int numSpecies_;

   real RectangleX1_;
   real RectangleY1_;
   real RectangleX2_;
   real RectangleY2_;
   
   real CircleX_;
   real CircleY_;
   real CircleR_;
};

inline real StaggeredGrid::p(const int x, const int y, Direction dir) {
  if (dir == NORTH) {
    if (isFluid(x,y+1)) {
      return p_(x,y+1);
    } else {
      return p_(x,y);
    }
  } else if (dir == SOUTH) {
    if (isFluid(x,y-1)) {
      return p_(x,y-1);
    } else {
      return p_(x,y);
    }
  } else if (dir == EAST) {
    if (isFluid(x+1,y)) {
      return p_(x+1,y);
    } else {
      return p_(x,y);
    }
  } else if (dir == WEST) { 
    if (isFluid(x-1,y)) {
      return p_(x-1,y);
    } else {
      return p_(x,y);
    }
  }
  CHECK(false)
  return 1.0/0.0;
}

inline real StaggeredGrid::u(const int x, const int y, Direction dir) {
  if (dir == NORTH) {
    if (!isFluid(x,y+1) != !isFluid(x+1,y+1)) {
      return 0.0;
    } else if (!isFluid(x,y+1) && !isFluid(x+1,y+1)) {
      return -u_(x,y);
    } else {
      return u_(x,y+1);
    }
  } else if (dir == SOUTH) {
    if (!isFluid(x,y-1) != !isFluid(x+1,y-1)) {
      return 0.0;
    } else if (!isFluid(x,y-1) && !isFluid(x+1,y-1)) {
      return -u_(x,y);
    } else {
      return u_(x,y-1);
    }
  } else if (dir == EAST) {
    if (isFluid(x+2,y)) {
      return u_(x+1,y);
    } else {
      return 0.0;
    }
  } else if (dir == WEST) { 
    if (isFluid(x-1,y)) {
      return u_(x-1,y);
    } else {
      return 0.0;
    }
  }
  CHECK(false)
  return 1.0/0.0;
}

inline real StaggeredGrid::v(const int x, const int y, Direction dir) {
  if (dir == NORTH) {
    if (isFluid(x,y+2)) {
      return v_(x,y+1);
    } else {
      return 0.0;
    }
  } else if (dir == SOUTH) {
    if (isFluid(x,y-1)) {
      return v_(x,y-1);
    } else {
      return 0.0;
    }
  } else if (dir == EAST) {
    if (!isFluid(x+1,y) != !isFluid(x+1,y+1)) {
      return 0.0;
    } else if (!isFluid(x+1,y) && !isFluid(x+1,y+1)) {
      return -v_(x,y);
    } else {
      return v_(x+1,y);
    }
  } else if (dir == WEST) { 
    if (!isFluid(x-1,y) != !isFluid(x-1,y+1)) {
      return 0.0;
    } else if (!isFluid(x-1,y) && !isFluid(x-1,y+1)) {
      return -v_(x,y);
    } else {
      return v_(x-1,y);
    }
  }
  CHECK(false)
  return 1.0/0.0;
}

inline bool StaggeredGrid::isFluid(const int x, const int y) {
   return !obstacleflags_(x,y);
}

inline void StaggeredGrid::setCellToObstacle(const int x, const int y) {
   obstacleflags_(x,y) = true;
   return;
}

inline int StaggeredGrid::getNumFluid() {
   int result = 0;
   for (int j = 1; j <= jmax_; j++) {   
     for (int i = 1; i <= imax_; i++) {
       isFluid(i,j) ? ++result : 0;              
     }
   }
   return result;
}

#endif //STAGGERED_GRID_HH

