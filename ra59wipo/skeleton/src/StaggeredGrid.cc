
#include <StaggeredGrid.hh>
#include <Debug.hh>
#include <iostream>


//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

StaggeredGrid::StaggeredGrid(int iimax, int jjmax, real ddx, real ddy) : p_(Array<real>(iimax+2, jjmax+2)),
    rhs_(Array<real>(iimax+2, jjmax+2)), u_(Array<real>(iimax+1, jjmax+2)), v_(Array<real>(iimax+2, jjmax+1)), f_(Array<real>(iimax+1,jjmax)), g_(Array<real>(iimax, jjmax+1)), obstacleflags_(Array<bool>(iimax+2, jjmax+2)), dx_(ddx), dy_(ddy), imax_(iimax), jmax_(jjmax) {
    obstacleflags_.fill(false);
}

StaggeredGrid::StaggeredGrid(const FileReader & configuration) {
    imax_ = configuration.getIntParameter("imax");
    jmax_ = configuration.getIntParameter("jmax");

    CHECK_MSG((imax_ >= 3) && (jmax_ >= 3), "Domain discretization too small. Check imax and jmax in your configuration.");
    
    real xLength = configuration.getRealParameter("xlength");
    real yLength = configuration.getRealParameter("ylength");
    CHECK_MSG((xLength > 0) && (yLength > 0), "Domain length too small. Check xlength and ylength in your configuration.");
    
    dx_ = xLength/imax_;
    dy_ = yLength/jmax_;

    p_ = Array<real>(imax_+2, jmax_+2);
    rhs_ = Array<real>(imax_+2, jmax_+2);
    u_ = Array<real>(imax_+1, jmax_+2);
    v_ = Array<real>(imax_+2, jmax_+1);
    f_ = Array<real>(imax_+1, jmax_+2);
    g_ = Array<real>(imax_+2, jmax_+1);
    obstacleflags_ = Array<bool>(imax_+2, jmax_+2);
    obstacleflags_.fill(false);

    numSpecies_ = configuration.getIntParameter("numSpecies");
    

    for(int i = 0; i < numSpecies_; ++i) {
        c_.push_back(Array<real>(imax_+2, jmax_+2));
        std::stringstream ss;
        ss << (i+1);
        lambda_.push_back(configuration.getRealParameter("lambda" + ss.str()));
        c(i).fill(configuration.getRealParameter("C" + ss.str() + "_INIT"));
    }
   
    RectangleX1_ = configuration.getRealParameter("RectangleX1");
    RectangleY1_ = configuration.getRealParameter("RectangleY1");
    RectangleX2_ = configuration.getRealParameter("RectangleX2");
    RectangleY2_ = configuration.getRealParameter("RectangleY2");

    CircleX_ = configuration.getRealParameter("CircleX");
    CircleY_ = configuration.getRealParameter("CircleY");
    CircleR_ = configuration.getRealParameter("CircleR");

    createRectangle(static_cast<int>(std::round(RectangleX1_/dx_)),
                    static_cast<int>(std::round(RectangleY1_/dy_)),
                    static_cast<int>(std::round(RectangleX2_/dx_)),
                    static_cast<int>(std::round(RectangleY2_/dy_)));

    createCircle(static_cast<int>(std::round(CircleX_/dx_)),
                 static_cast<int>(std::round(CircleY_/dy_)),
                 static_cast<int>(std::round(CircleR_/std::max(dx_,dy_))));
    
    saveObstacleFieldToPNG("obstaclefield.png");
}

void StaggeredGrid::registerModule(FileReader& configuration) {
    configuration.registerRealParameter("RectangleX1",-2.0);
    configuration.registerRealParameter("RectangleY1",-2.0);
    configuration.registerRealParameter("RectangleX2",-1.0);
    configuration.registerRealParameter("RectangleY2",-1.0);
    configuration.registerRealParameter("CircleX",0.0);
    configuration.registerRealParameter("CircleY",0.0);
    configuration.registerRealParameter("CircleR",-1.0);
    return;
}

void StaggeredGrid::createRectangle(int x1, int y1, int x2, int y2) {
   if (std::min(std::min(std::min(x1,y1),x2),y2) < 0) return;
   x1 = std::min(std::max(x1,1),imax_+1);
   x2 = std::min(std::max(x2,1),imax_+1);  
   y1 = std::min(std::max(y1,1),jmax_+1);
   y2 = std::min(std::max(y2,1),jmax_+1);
   for (int j = y1; j <= y2; j++) {
     for (int i = x1; i <= x2; i++) {
       setCellToObstacle(i,j);
     }
   }
   return;
}

void StaggeredGrid::createCircle(int x, int y, int r) {
   if (r < 0) return;
   for (int j = 1; j <= jmax_; j++) {
     for (int i = 1; i <= imax_; i++) {
       if ((i-x) * (i-x) + (j-y) * (j-y) <= r * r) {
         setCellToObstacle(i,j);
       }
     }
   }
   return;
}

void StaggeredGrid::saveObstacleFieldToPNG(std::string filename) {
   static const unsigned char maxVal = std::numeric_limits<unsigned char>::max();
   GrayScaleImage png(imax_+2, jmax_+2);
   for (int j = 0; j <= jmax_+1; j++) {
    for (int i = 0; i <= imax_+1; i++) {
      png(i,j) = isFluid(i,j) ? maxVal : 0;
    }
   }
   png.save(filename);
   return; 
}

void StaggeredGrid::loadObstacleFieldFromPNG(std::string filename) {
   GrayScaleImage png(filename);
   CHECK_MSG((png.width() == imax_+2) && (png.height() == jmax_+2), "Image dimension does not fit the simulation domain.")
   for (int j = 1; j <= jmax_; j++) {
     for (int i = 1; i <= imax_; i++) {
       if (png(i,j) == 0.0) {
         setCellToObstacle(i,j);
       }
     }
   }
   return;
}
