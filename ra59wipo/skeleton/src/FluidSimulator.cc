#include "FluidSimulator.hh"
#include <vector>
#include <iostream>
using namespace std;
//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

void FluidSimulator::registerModule(FileReader & params) {
   params.registerStringParameter("boundary_condition_S","NOSLIP");
   params.registerStringParameter("boundary_condition_N","NOSLIP");
   params.registerStringParameter("boundary_condition_W","NOSLIP");
   params.registerStringParameter("boundary_condition_E","NOSLIP");
   params.registerRealParameter("boundary_velocity_S",0.0);
   params.registerRealParameter("boundary_velocity_N",0.0);
   params.registerRealParameter("boundary_velocity_W",0.0);
   params.registerRealParameter("boundary_velocity_E",0.0);
   params.registerRealParameter("safetyfactor",0.9);
   params.registerRealParameter("U_INIT",0.0);
   params.registerRealParameter("V_INIT",0.0);
   params.registerRealParameter("P_INIT",1.0);
   params.registerIntParameter("outputinterval",1);
   StaggeredGrid::registerModule(params);
   SORSolver::registerModule(params);
}

FluidSimulator::FluidSimulator(const FileReader & conf) : sg_(StaggeredGrid(conf)), sor_(SORSolver(conf)), vtkWriter_(sg_, conf.getStringParameter("name"), true, true) {
    dt_ = conf.getRealParameter("dt");
    CHECK_MSG((dt_ >= 0), "dt is negative. Check your configuration.");
    Re_ = conf.getRealParameter("Re");
    CHECK_MSG((Re_ >= 0), "Re is negative. Check your configuration.");
    gamma_ = conf.getRealParameter("gamma");
    CHECK_MSG(((gamma_ >= 0) && (gamma_ <=1)), "gamma is not in [0,1]. Check your configuration.");
    gx_ = conf.getRealParameter("gx");
    gy_ = conf.getRealParameter("gy");
    safetyfactor_ = conf.getRealParameter("safetyfactor");
    boundary_condition_S_ = conf.getStringParameter("boundary_condition_S");
    boundary_condition_N_ = conf.getStringParameter("boundary_condition_N");
    boundary_condition_W_ = conf.getStringParameter("boundary_condition_W");
    boundary_condition_E_ = conf.getStringParameter("boundary_condition_E");
    boundary_velocity_S_ = conf.getRealParameter("boundary_velocity_S");
    boundary_velocity_N_ = conf.getRealParameter("boundary_velocity_N");
    boundary_velocity_W_ = conf.getRealParameter("boundary_velocity_W");
    boundary_velocity_E_ = conf.getRealParameter("boundary_velocity_E");
    U_INIT_ = conf.getRealParameter("U_INIT");
    V_INIT_ = conf.getRealParameter("V_INIT");
    P_INIT_ = conf.getRealParameter("P_INIT");

    outputinterval_ = conf.getIntParameter("outputinterval");
}

void FluidSimulator::determineNextDT ( real const & limit ) {
    if (safetyfactor_ > 0) {
        Array<real> & u = sg_.u();
        Array<real> & v = sg_.v();
        auto minimax_u = u.minmax();
        auto minimax_v = v.minmax();
        auto dx = sg_.dx();
        auto dy = sg_.dy();
        dt_ = safetyfactor_ * std::min(std::min(Re_/(2.0 * (1.0/(dx * dx) + 1.0/(dy * dy))),
                                       dx / std::max(fabs(minimax_u.first), minimax_u.second)),
                                       dy / std::max(fabs(minimax_v.first), minimax_v.second));
    }
    return;
}

void FluidSimulator::updateVelocities() {
    Array<real> & u = sg_.u();
    Array<real> & v = sg_.v();
    Array<real> & f = sg_.f();
    Array<real> & g = sg_.g();
    Array<real> & p = sg_.p();
    real invDx = 1.0/sg_.dx();
    real invDy = 1.0/sg_.dy();
    int imax = sg_.imax();
    int jmax = sg_.jmax();
    
    for (int j = 1; j <= jmax; j++) {
        for (int i = 1; i <= imax - 1; i++) {
            if (!sg_.isFluid(i,j)) continue;        
            u(i,j) = f(i,j) - dt_ * invDx * (sg_.p(i,j,EAST) - p(i,j));
        }
    }
    
    for (int j = 1; j <= jmax - 1; j++) {
        for (int i = 1; i <= imax; i++) {
            if (!sg_.isFluid(i,j)) continue;        
            v(i,j) = g(i,j) - dt_ * invDy * (sg_.p(i,j,NORTH) - p(i,j));
        }
    }
    
    return;    
}

void FluidSimulator::composeRHS() {
    Array<real> & f = sg_.f();
    Array<real> & g = sg_.g();
    Array<real> & rhs = sg_.rhs();
    
    int imax = sg_.imax();
    int jmax = sg_.jmax();
    
    real invDt = 1.0/dt_;
    real invDx = 1.0/sg_.dx();
    real invDy = 1.0/sg_.dy();
    
    real average = 0.0;
    for (int j = 1; j <= jmax; j++) {
        for (int i = 1; i <= imax; i++) {
            if (!sg_.isFluid(i,j)) continue;
            rhs(i,j) = invDt * ( invDx * (f(i,j) - f(i-1,j)) + invDy * (g(i,j) - g(i,j-1)) );
            average += rhs(i,j);
        }
    }
    
    average /= sg_.getNumFluid();
    for (int j = 1; j <= jmax; j++) {
        for (int i = 1; i <= imax; i++) {
            if (!sg_.isFluid(i,j)) continue;
            rhs(i,j) -= average;
        }
    }
    
    return;
}

void FluidSimulator::refreshBoundaries() {
    Array<real> & u = sg_.u();
    Array<real> & v = sg_.v();
    
    int imax = sg_.imax();
    int jmax = sg_.jmax();
    if (boundary_condition_S_.compare("NOSLIP") == 0) {
        for (int i = 1; i <= imax; i++) {
            v(i,0) = 0.0;
        }
        for (int i = 1; i <= imax; i++) {
            u(i,0) = 2.0 * boundary_velocity_S_ - u(i,1);
        }
    } else if (boundary_condition_S_.compare("INFLOW") == 0) {
        for (int i = 1; i <= imax; i++) {
            v(i,0) = boundary_velocity_S_;
        }
        for (int i = 1; i <= imax; i++) {
            u(i,0) = 0.0;
        }
    } else if (boundary_condition_S_.compare("OUTFLOW") == 0) {
        CHECK_MSG((boundary_velocity_S_ == 0.0), "OUTFLOW boundary velocity must not be set.");
        for (int i = 1; i <= imax; i++) {
            v(i,0) = v(i,1);
        }
        for (int i = 1; i <= imax; i++) {
            u(i,0) = u(i,1);
        }	
    } else {
        CHECK_MSG(false, boundary_condition_S_);
    }
    
    if (boundary_condition_N_.compare("NOSLIP") == 0) {
        for (int i = 1; i <= imax; i++) {
            v(i,jmax) = 0.0;
        }
        for (int i = 1; i <= imax; i++) {
            u(i,jmax+1) = 2.0 * boundary_velocity_N_ - u(i,jmax);
        }
    } else if (boundary_condition_N_.compare("INFLOW") == 0) {
        for (int i = 1; i <= imax; i++) {
            v(i,jmax) = boundary_velocity_N_;
        }
        for (int i = 1; i <= imax; i++) {
            u(i,jmax+1) = 0.0;
        }
    } else if (boundary_condition_N_.compare("OUTFLOW") == 0) {
        CHECK_MSG((boundary_velocity_N_ == 0.0), "OUTFLOW boundary velocity must not be set.");
	for (int i = 1; i <= imax; i++) {
            v(i,jmax) = v(i,jmax-1);
        }
        for (int i = 1; i <= imax; i++) {
            u(i,jmax+1) = u(i,jmax);
        }
    } else {
        CHECK_MSG(false, boundary_condition_N_);
    }
    
    if (boundary_condition_E_.compare("NOSLIP") == 0) {
        for (int j = 1; j <= jmax; j++) {
            u(imax,j) = 0.0;
        }
        for (int j = 1; j <= jmax; j++) {
            v(imax+1,j) = 2.0 * boundary_velocity_E_ - v(imax,j);
        }
        
    } else if (boundary_condition_E_.compare("INFLOW") == 0) {
        for (int j = 1; j <= jmax; j++) {
            u(imax,j) = boundary_velocity_E_;
        }
        for (int j = 1; j <= jmax; j++) {
            v(imax+1,j) = 0.0;
        }
    } else if (boundary_condition_E_.compare("OUTFLOW") == 0) {
        CHECK_MSG((boundary_velocity_E_ == 0.0), "OUTFLOW boundary velocity must not be set.");
	for (int j = 1; j <= jmax; j++) {
            u(imax,j) = u(imax-1,j);
        }
        for (int j = 1; j <= jmax; j++) {
            v(imax+1,j) = v(imax,j);
        }
    } else {
        CHECK_MSG(false, boundary_condition_E_);
    }

    if (boundary_condition_W_.compare("NOSLIP") == 0) {
        for (int j = 1; j <= jmax; j++) {
            u(0,j) = 0.0;
        }
        for (int j = 1; j <= jmax; j++) {
            v(0,j) = 2.0 * boundary_velocity_W_ - v(1,j);
        }
    } else if (boundary_condition_W_.compare("INFLOW") == 0) {
        for (int j = 1; j <= jmax; j++) {
            u(0,j) = boundary_velocity_W_;
        }
        for (int j = 1; j <= jmax; j++) {
            v(0,j) = 0.0;
        }
    } else if (boundary_condition_W_.compare("OUTFLOW") == 0) {
        CHECK_MSG((boundary_velocity_W_ != 0.0), "OUTFLOW boundary velocity must not be set.");
	for (int j = 1; j <= jmax; j++) {
            u(0,j) = u(1,j);
        }
        for (int j = 1; j <= jmax; j++) {
            v(0,j) = v(1,j);
        }
    } else {
        CHECK_MSG(false, boundary_condition_W_);
    }
    
    return;
}

void FluidSimulator::updateCBoundaries(){
  
}

void FluidSimulator::updateC() {

    real invDx = 1.0/sg_.dx();
    real invDy = 1.0/sg_.dy();

    real invDx2 = invDx*invDx;
    real invDy2 = invDy*invDy;
    //real invDt = 1.0;

    Array<real> & u = sg_.u();
    Array<real> & v = sg_.v();

    for(int k = 0; k < sg_.numSpecies(); ++k){

        Array<real> & c = sg_.c(k);
        
        for(int j = 1; j <= sg_.jmax(); ++j) {
            for(int i = 1; i <= sg_.imax(); ++i) { 
                if (!sg_.isFluid(i,j)) continue;  
      
                real ducdx = invDx * 0.5 * (u(i,j)*(c(i,j)+sg_.c(i,j,EAST,k)) - sg_.u(i,j,WEST)*(sg_.c(i,j,WEST,k)+c(i,j))) + gamma_*invDx * 0.5 * (fabs(u(i,j))*(c(i,j)-sg_.c(i,j,EAST,k)) - fabs(sg_.u(i,j,WEST))*(sg_.c(i,j,WEST,k)-c(i,j)));

                real dvcdy = invDy * 0.5 * (v(i,j)*(c(i,j)+sg_.c(i,j,NORTH,k)) - sg_.v(i,j,SOUTH)*(sg_.c(i,j,SOUTH,k)+c(i,j))) + gamma_*invDy * 0.5 * (fabs(v(i,j))*(c(i,j)-sg_.c(i,j,NORTH,k)) - fabs(sg_.v(i,j,SOUTH))*(sg_.c(i,j,SOUTH,k)-c(i,j)));
           
                real d2cdx2 = invDx2 * (sg_.c(i,j,EAST,k) - 2.0*c(i,j) + sg_.c(i,j,WEST,k));

                real d2cdy2 = invDy2 * (sg_.c(i,j,NORTH,k) - 2.0*c(i,j) + sg_.c(i,j,SOUTH,k));

                c(i,j) = c(i,j) + dt_ * (sg_.lambda(k)*(d2cdx2+d2cdy2)*c(i,j) /*TODO: + q() */ - ducdx - dvcdy);

            }
        }
    }


}

void FluidSimulator::computeFG() {
    Array<real> & f = sg_.f();
    Array<real> & g = sg_.g();
    Array<real> & u = sg_.u();
    Array<real> & v = sg_.v();
    real reynoldsFactor = 1.0/Re_;
    real dx = sg_.dx();
    real dxFactor = 1.0/(dx * dx);
    real invDx = 1.0/dx;
    real dy = sg_.dy();
    real dyFactor = 1.0/(dy * dy);
    real invDy = 1.0/dy;

    int imax = sg_.imax();
    int jmax = sg_.jmax();

    // Set Boundary Conditions
    for (int j = 1; j <= jmax; j++) {
        f(0,j) = u(0,j);
        f(imax,j) = u(imax,j);
    }
    
    for (int i = 1; i <= imax; i++) {
        g(i,0) = v(i,0);
        g(i,jmax) = v(i,jmax);
    }
    
    // Iterate over the domain
    for (int j = 1; j <= jmax; j++) {
        for (int i = 1; i < imax; i++) {
            if (!sg_.isFluid(i,j) || !sg_.isFluid(i,j+1)) continue;        
            real ReLaplaceU = reynoldsFactor * (dxFactor * (sg_.u(i,j,EAST) - 2.0 * u(i,j) + sg_.u(i,j,WEST)) +
                    dyFactor * (sg_.u(i,j,NORTH) - 2.0 * u(i,j) + sg_.u(i,j,SOUTH)));
            real duudx = invDx * 0.25 * ( ((u(i,j)+sg_.u(i,j,EAST)) * (u(i,j)+sg_.u(i,j,EAST)) - (sg_.u(i,j,WEST)+u(i,j)) * (sg_.u(i,j,WEST)+u(i,j)) )
                    + gamma_ * ( fabs(u(i,j)+sg_.u(i,j,EAST)) * (u(i,j)-sg_.u(i,j,EAST)) - fabs(sg_.u(i,j,WEST) + u(i,j)) * (sg_.u(i,j,WEST) - u(i,j)) ) );
            real duvdy = invDy * 0.25 * ( ( (sg_.v(i,j-1,NORTH) + sg_.v(i+1,j-1,NORTH)) * (u(i,j) + sg_.u(i,j,NORTH)) - (sg_.v(i,j,SOUTH) + sg_.v(i+1,j,SOUTH)) * (sg_.u(i,j,SOUTH) + u(i,j)))
                    + gamma_ * ( fabs(sg_.v(i,j-1,NORTH) + sg_.v(i+1,j-1,NORTH)) * (u(i,j) - sg_.u(i,j,NORTH)) - fabs(sg_.v(i,j,SOUTH) + sg_.v(i+1,j,SOUTH)) * (sg_.u(i,j,SOUTH) - u(i,j))) );
            f(i,j) = u(i,j) + dt_ * (ReLaplaceU - duudx - duvdy + gx_);
    
        }
    }
        
    for (int j = 1; j < jmax; j++) {
        for (int i = 1; i <= imax; i++) {
            if (!sg_.isFluid(i,j) || !sg_.isFluid(i+1,j)) continue;        
            real ReLaplaceV = reynoldsFactor * (dxFactor * (sg_.v(i,j,EAST) - 2.0 * v(i,j) + (sg_.v(i,j,WEST))) +
                    dyFactor * (sg_.v(i,j,NORTH) - 2.0 * v(i,j) + sg_.v(i,j,SOUTH)));
            real duvdx = 0.25 * invDx * ( ( (sg_.u(i-1,j,EAST) + sg_.u(i-1,j+1,EAST)) * (v(i,j) + sg_.v(i,j,EAST)) - (sg_.u(i,j,WEST) + sg_.u(i,j+1,WEST)) * (sg_.v(i,j,WEST) + v(i,j)) )
                    + gamma_ * ( fabs(u(i,j) + sg_.u(i-1,j+1,EAST)) *(v(i,j) - sg_.v(i,j,EAST)) - fabs(sg_.u(i,j,WEST) + sg_.u(i,j+1,WEST)) * (sg_.v(i,j,WEST) - v(i,j)) ));
            real dvvdy = 0.25 * invDy * ( ( (v(i,j)+sg_.v(i,j,NORTH))*(v(i,j)+sg_.v(i,j,NORTH)) - (sg_.v(i,j,SOUTH)+v(i,j))*(sg_.v(i,j,SOUTH)+v(i,j)) ) -
                    + gamma_ * ( fabs(v(i,j)+sg_.v(i,j,NORTH))*(v(i,j)-sg_.v(i,j,NORTH)) - fabs(sg_.v(i,j,SOUTH)+v(i,j))*(sg_.v(i,j,SOUTH)-v(i,j)) ));

            g(i,j) = v(i,j) + dt_ * (ReLaplaceV - duvdx - dvvdy + gy_);
        }
    }
}

void FluidSimulator::assignInitialValues() {
    sg_.u().fill(U_INIT_);
    sg_.v().fill(V_INIT_);
    sg_.p().fill(P_INIT_);
    sg_.f().fill(0.0);
    sg_.g().fill(0.0);
    sg_.rhs().fill(0.0);    
    return;
}

void FluidSimulator::simulateGeneral(std::function<bool(real,unsigned int)> criterion) {
    real t = 0.0;
    unsigned int n = 0;
    assignInitialValues();
    
    while (criterion(t,n)) {
        {
          std::stringstream is;
          is  << "Simulation of timestep: " << n << " time: " << t;
          PRG_MSG(is.str())
        }
        determineNextDT(0.0);
        refreshBoundaries();
        computeFG();
        composeRHS();
        sor_.solve(sg_);
        //CHECK(sor_.solve(sg_));
        updateVelocities();
        {
          if ((outputinterval_ != 0) && (n % outputinterval_ == 0)) {  
            vtkWriter_.write();
            std::stringstream is;
            is  << "Writing VTK file.";
            PRG_MSG(is.str())
          }
        }
        t += dt_;
        n += 1;
    }    
}

void FluidSimulator::simulate(real duration) {
    simulateGeneral([duration](real t, int n) -> bool {return t <= duration;});
}

void FluidSimulator::simulateTimeStepCount(unsigned int nrOfTimeSteps) {
    simulateGeneral([nrOfTimeSteps](real t, unsigned int n) -> bool {return n <= nrOfTimeSteps;});
}
