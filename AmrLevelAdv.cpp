
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCUtil.H>
#include <iostream>
#include <fstream>
#include <string>


using namespace amrex;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9; // Default value - can be overwritten in settings file
int      AmrLevelAdv::do_reflux       = 1;  

int      AmrLevelAdv::NUM_STATE       = 4;  // Changed to 4 variables for the 2D Euler equations - mc2246 
int      AmrLevelAdv::NUM_GROW        = 2;  // Changed to 2 ghost cells for transmissive boundary conditions - mc2246

// Added function to compute euler flux in x-direction - mc2246
// Euler flux in x-direction
std::array<double,4> feuler_x(Array4<Real> const &arr, int i, int j, int k){
	
	const double gamma = 1.4;
	std::array<double, 4> flux_x;
	
	double rho = arr(i,j,k,0);
	double momx = arr(i,j,k,1);
	double momy = arr(i,j,k,2);
	double E = arr(i,j,k,3);
	
	double vx = arr(i,j,k,1)/arr(i,j,k,0);
	double vy = arr(i,j,k,2)/arr(i,j,k,0);
	double p = (gamma - 1.0)*(E - 0.5*rho*(vx*vx + vy*vy));
	
	flux_x[0] = momx;
	flux_x[1] = momx*vx + p; 
	flux_x[2] = momx*vy;
	flux_x[3] = (E + p)*vx;
	
	return flux_x;
}

// Added function to compute euler flux in y-direction - mc2246
// Euler flux in y-direction
std::array<double,4> feuler_y(Array4<Real> const &arr, int i, int j, int k){
	
	const double gamma = 1.4;
	std::array<double, 4> flux_y;
	
	double rho = arr(i,j,k,0);
	double momx = arr(i,j,k,1);
	double momy = arr(i,j,k,2);
	double E = arr(i,j,k,3);
	
	double vx = arr(i,j,k,1)/arr(i,j,k,0);
	double vy = arr(i,j,k,2)/arr(i,j,k,0);
	double p = (gamma - 1.0)*(E - 0.5*rho*(vx*vx + vy*vy));
	
	flux_y[0] = momy;
	flux_y[1] = momy*vx; //check this 
	flux_y[2] = momy*vy + p;
	flux_y[3] = (E + p)*vy;
	
	return flux_y;
}

// Added function to use in MUSCL-Hancock Method - mc2246
// Compute Delta in x-direction
double compute_deltax(Array4<Real> const &arr, const int i, const int j, const int k, const int n, const double omega){
	
	double delta;
	delta =  0.5*(1.0 + omega)*(arr(i,j,k,n) - arr(i-1,j,k,n)) + 0.5*(1.0 - omega)*(arr(i+1,j,k,n) - arr(i,j,k,n)); 
	
	return delta;
}

// Added function to use in MUSCL-Hancock Method - mc2246
// Compute Delta in y-direction
double compute_deltay(Array4<Real> const &arr, const int i, const int j, const int k, const int n, const double omega){
	
	double delta;
	delta =  0.5*(1.0 + omega)*(arr(i,j,k,n) - arr(i,j-1,k,n)) + 0.5*(1.0 - omega)*(arr(i,j+1,k,n) - arr(i,j,k,n)); 
	
	return delta;
}

// Added function to use in MUSCL-Hancock Method - mc2246
// Compute quantity r in x-direction
double compute_rx(Array4<Real> const &arr, const int i, const int j, const int k, const int n){
	
	double r = 0.0;
	double r_num = arr(i,j,k,n) - arr(i-1,j,k,n);
	double r_denom = arr(i+1,j,k,n) - arr(i,j,k,n);
	 
	if (r_denom == 0){
		r = 0.0;
	}
	else{
		r = r_num/r_denom;
	}
	
	return r;
}

// Added function to use in MUSCL-Hancock Method - mc2246
// Compute quantity r in y-direction
double compute_ry(Array4<Real> const &arr, const int i, const int j, const int k, const int n){
	
	double r = 0.0;
	double r_num = arr(i,j,k,n) - arr(i,j-1,k,n);
	double r_denom = arr(i,j+1,k,n) - arr(i,j,k,n);
	 
	if (r_denom == 0){
		r = 0.0;
	}
	else{
		r = r_num/r_denom;
	}
	
	return r;
}

// Added slope limiter function to use in MUSCL-Hancock Method - mc2246
// Compute Minbee Slope limiter
double Minbeelimiter(const double &r){
	
	double minbee_psi = 0.0;
	
	if (r <= 0){
		minbee_psi = 0.0;		
	}
	else if (0<r && r<=1){
		minbee_psi = r;
	}
	else if (r>1){
		double epsilonR = 2.0/(1.0 + r);
		minbee_psi = std::min(1.0, epsilonR); 
	}
	
	return minbee_psi;
}

// Added HLL flux function to check code is working properly before implementing HLLC - mc2246
// HLL flux in x-direction
std::array<double,4> getHLLflux_x(Array4<Real> const &UR_nph_x, Array4<Real> const &UL_nph_x, int i, int j, int k){
	//Left state = UR_nph_x(i-1,j,k)
	//Right state = UL_nph_x(i,j,k)
	
	const double gamma = 1.4;
	
	// Left state	
	Real rho_L = UR_nph_x(i-1,j,k,0);
	//Real momx_L = UR_nph_x(i-1,j,k,1);
	//Real momy_L = UR_nph_x(i-1,j,k,2);
	Real E_L = UR_nph_x(i-1,j,k,3);
	Real vx_L = UR_nph_x(i-1,j,k,1)/UR_nph_x(i-1,j,k,0);
	Real vy_L = UR_nph_x(i-1,j,k,2)/UR_nph_x(i-1,j,k,0);
	Real p_L = (gamma - 1.0)*(E_L - 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L));			
	Real Cs_L = sqrt(gamma*(p_L/rho_L));
	
	//Right state 
	Real rho_R = UL_nph_x(i,j,k,0);
	//Real momx_R = UL_nph_x(i,j,k,1);
	//Real momy_R = UL_nph_x(i,j,k,2);
	Real E_R = UL_nph_x(i,j,k,3);
	Real vx_R = UL_nph_x(i,j,k,1)/UL_nph_x(i,j,k,0);
	Real vy_R = UL_nph_x(i,j,k,2)/UL_nph_x(i,j,k,0);
	Real p_R = (gamma - 1.0)*(E_R - 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R));			
	Real Cs_R = sqrt(gamma*(p_R/rho_R));
	
	//Wave speed estimates //might have to change these
	Real SL = std::min(vx_L - Cs_L, vx_R - Cs_R);
	Real SR = std::max(vx_L + Cs_L, vx_R + Cs_R);
	
	//Get euler flux
	std::array<double,4> fx_L = feuler_x(UR_nph_x, i-1, j, k);
	std::array<double,4> fx_R = feuler_x(UL_nph_x, i, j, k);
	
	//Compute HLL flux 
	std::array<double, 4> HLLflux_x; 
	for (int n = 0; n < 4; n++){
		HLLflux_x[n] = (SR*fx_L[n] - SL*fx_R[n] + SL*SR*(UL_nph_x(i,j,k,n) - UR_nph_x(i-1,j,k,n)))/(SR - SL);
	}
	
	std::array<double, 4> finalflux;
	if (0 <= SL){
		finalflux = fx_L;
	}
	else if(SL < 0 && 0< SR){
		finalflux = HLLflux_x;
	}
	else if(SR <= 0){
		finalflux = fx_R; 
	}
		
	return finalflux;
}

// Added HLL flux function to check code is working properly before implementing HLLC - mc2246
// HLL flux in y-direction
std::array<double,4> getHLLflux_y(Array4<Real> const &UR_nph_y, Array4<Real> const &UL_nph_y, int i, int j, int k){
	//Left state = UR_nph_y(i,j-1,k)
	//Right state = UL_nph_y(i,j,k)
	
	const double gamma = 1.4;
	
	// Left state	
	Real rho_L = UR_nph_y(i,j-1,k,0);
	//Real momx_L = UR_nph_x(i,j-1,k,1);
	//Real momy_L = UR_nph_x(i,j-1,k,2);
	Real E_L = UR_nph_y(i,j-1,k,3);
	Real vx_L = UR_nph_y(i,j-1,k,1)/UR_nph_y(i,j-1,k,0);
	Real vy_L = UR_nph_y(i,j-1,k,2)/UR_nph_y(i,j-1,k,0);
	Real p_L = (gamma - 1.0)*(E_L - 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L));			
	Real Cs_L = sqrt(gamma*(p_L/rho_L));
	
	//Right state 
	Real rho_R = UL_nph_y(i,j,k,0);
	//Real momx_R = UL_nph_x(i+2,j,k,1);
	//Real momy_R = UL_nph_x(i+2,j,k,2);
	Real E_R = UL_nph_y(i,j,k,3);
	Real vx_R = UL_nph_y(i,j,k,1)/UL_nph_y(i,j,k,0);
	Real vy_R = UL_nph_y(i,j,k,2)/UL_nph_y(i,j,k,0);
	Real p_R = (gamma - 1.0)*(E_R - 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R));			
	Real Cs_R = sqrt(gamma*(p_R/rho_R));
	
	//Wave speed estimates //might have to change these
	Real SL = std::min(vy_L - Cs_L, vy_R - Cs_R);
	Real SR = std::max(vy_L + Cs_L, vy_R + Cs_R);
	
	//Get euler flux
	std::array<double,4> fy_L = feuler_y(UR_nph_y, i, j-1, k);
	std::array<double,4> fy_R = feuler_y(UL_nph_y, i, j, k);
	
	//Compute HLL flux 
	std::array<double, 4> HLLflux_y; 
	for (int n = 0; n < 4; n++){
		HLLflux_y[n] = (SR*fy_L[n] - SL*fy_R[n] + SL*SR*(UL_nph_y(i,j,k,n) - UR_nph_y(i,j-1,k,n)))/(SR - SL);
	}
	
	std::array<double, 4> finalflux;
	if (0 <= SL){
		finalflux = fy_L;
	}
	else if(SL < 0 && 0< SR){
		finalflux = HLLflux_y;
	}
	else if(SR <= 0){
		finalflux = fy_R; 
	}
		
	return finalflux;
}

// Added HLLC flux function - mc2246
// HLLC flux in x-direction
std::array<double,4> getHLLCflux_x(Array4<Real> const &UR_nph_x, Array4<Real> const &UL_nph_x, int i, int j, int k){
	//Left state = UR_nph_x(i-1,j,k)
	//Right state = UL_nph_x(i,j,k)
	
	const double gamma = 1.4;
	
	// Left state	
	Real rho_L = UR_nph_x(i-1,j,k,0);
	//Real momx_L = UR_nph_x(i-1,j,k,1);
	//Real momy_L = UR_nph_x(i-1,j,k,2);
	Real E_L = UR_nph_x(i-1,j,k,3);
	Real vx_L = UR_nph_x(i-1,j,k,1)/UR_nph_x(i-1,j,k,0);
	Real vy_L = UR_nph_x(i-1,j,k,2)/UR_nph_x(i-1,j,k,0);
	Real p_L = (gamma - 1.0)*(E_L - 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L));			
	Real Cs_L = sqrt(gamma*(p_L/rho_L));
	
	//Right state 
	Real rho_R = UL_nph_x(i,j,k,0);
	//Real momx_R = UL_nph_x(i,j,k,1);
	//Real momy_R = UL_nph_x(i,j,k,2);
	Real E_R = UL_nph_x(i,j,k,3);
	Real vx_R = UL_nph_x(i,j,k,1)/UL_nph_x(i,j,k,0);
	Real vy_R = UL_nph_x(i,j,k,2)/UL_nph_x(i,j,k,0);
	Real p_R = (gamma - 1.0)*(E_R - 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R));			
	Real Cs_R = sqrt(gamma*(p_R/rho_R));
	
	//Wave speed estimates //might have to change these
	
	Real Splus = std::max(fabs(vx_L)+Cs_L, fabs(vx_R)+Cs_R);
	Real SL = -Splus;
	Real SR = Splus;
	
	//Real SL = std::min(vx_L - Cs_L, vx_R - Cs_R);
	//Real SR = std::max(vx_L + Cs_L, vx_R + Cs_R);
	Real Sstar = (p_R - p_L + rho_L*vx_L*(SL - vx_L) - rho_R*vx_R*(SR - vx_R))/(rho_L*(SL - vx_L)-rho_R*(SR - vx_R));
	
	//Get euler flux
	std::array<double,4> fx_L = feuler_x(UR_nph_x, i-1, j, k);
	std::array<double,4> fx_R = feuler_x(UL_nph_x, i, j, k);
	
	//Compute HLLC left state
	std::array<double,4> UL_HLLC;
	UL_HLLC[0] = rho_L*((SL - vx_L)/(SL - Sstar));
	UL_HLLC[1] = UL_HLLC[0]*Sstar;
	UL_HLLC[2] = UL_HLLC[0]*vy_L;
	UL_HLLC[3] = UL_HLLC[0]*(E_L/rho_L + (Sstar - vx_L)*(Sstar + p_L/(rho_L*(SL - vx_L))));
	
	//Compute HLLC right state
	std::array<double,4> UR_HLLC;
	UR_HLLC[0] = rho_R*((SR - vx_R)/(SR - Sstar));
	UR_HLLC[1] = UR_HLLC[0]*Sstar;
	UR_HLLC[2] = UR_HLLC[0]*vy_R;
	UR_HLLC[3] = UR_HLLC[0]*(E_R/rho_R + (Sstar - vx_R)*(Sstar + p_R/(rho_R*(SR - vx_R))));
	
	//Compute HLL flux 
	std::array<double, 4> fL_HLLC_x;
	std::array<double, 4> fR_HLLC_x; 
	for (int n = 0; n < 4; n++){
		fL_HLLC_x[n] = fx_L[n] + SL*(UL_HLLC[n] - UR_nph_x(i-1,j,k,n));
		fR_HLLC_x[n] = fx_R[n] + SR*(UR_HLLC[n] - UL_nph_x(i,j,k,n));
	}
	
	std::array<double, 4> finalflux;
	if (0 <= SL){
		finalflux = fx_L;
	}
	else if(SL < 0 && 0 <= Sstar){
		finalflux = fL_HLLC_x;
	}
	else if(Sstar < 0 &&  0 <= SR){
		finalflux = fR_HLLC_x;
	}
	else if(SR < 0){
		finalflux = fx_R; 
	}
		
	return finalflux;
}

// Added HLLC flux function - mc2246
// HLLC flux in y-direction
std::array<double,4> getHLLCflux_y(Array4<Real> const &UR_nph_y, Array4<Real> const &UL_nph_y, int i, int j, int k){
	//Left state = UR_nph_y(i,j-1,k)
	//Right state = UL_nph_y(i,j,k)
	
	const double gamma = 1.4;
	
	// Left state	
	Real rho_L = UR_nph_y(i,j-1,k,0);
	//Real momx_L = UR_nph_x(i,j-1,k,1);
	//Real momy_L = UR_nph_x(i,j-1,k,2);
	Real E_L = UR_nph_y(i,j-1,k,3);
	Real vx_L = UR_nph_y(i,j-1,k,1)/UR_nph_y(i,j-1,k,0);
	Real vy_L = UR_nph_y(i,j-1,k,2)/UR_nph_y(i,j-1,k,0);
	Real p_L = (gamma - 1.0)*(E_L - 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L));			
	Real Cs_L = sqrt(gamma*(p_L/rho_L));
	
	//Right state 
	Real rho_R = UL_nph_y(i,j,k,0);
	//Real momx_R = UL_nph_x(i+2,j,k,1);
	//Real momy_R = UL_nph_x(i+2,j,k,2);
	Real E_R = UL_nph_y(i,j,k,3);
	Real vx_R = UL_nph_y(i,j,k,1)/UL_nph_y(i,j,k,0);
	Real vy_R = UL_nph_y(i,j,k,2)/UL_nph_y(i,j,k,0);
	Real p_R = (gamma - 1.0)*(E_R - 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R));			
	Real Cs_R = sqrt(gamma*(p_R/rho_R));
	
	//Wave speed estimates //might have to change these
	//Real SL = std::min(vy_L - Cs_L, vy_R - Cs_R);
	//Real SR = std::max(vy_L + Cs_L, vy_R + Cs_R);
	
	Real Splus = std::max(fabs(vy_L)+Cs_L, fabs(vy_R)+Cs_R);
	Real SL = -Splus;
	Real SR = Splus;
	
	Real Sstar = (p_R - p_L + rho_L*vy_L*(SL - vy_L) - rho_R*vy_R*(SR - vy_R))/(rho_L*(SL - vy_L)-rho_R*(SR - vy_R));
	
	//Get euler flux
	std::array<double,4> fy_L = feuler_y(UR_nph_y, i, j-1, k);
	std::array<double,4> fy_R = feuler_y(UL_nph_y, i, j, k);
	
	//Compute HLLC left state
	std::array<double,4> UL_HLLC;
	UL_HLLC[0] = rho_L*((SL - vy_L)/(SL - Sstar));
	UL_HLLC[1] = UL_HLLC[0]*vx_L;
	UL_HLLC[2] = UL_HLLC[0]*Sstar;
	UL_HLLC[3] = UL_HLLC[0]*(E_L/rho_L + (Sstar - vy_L)*(Sstar + p_L/(rho_L*(SL - vy_L))));
	
	//Compute HLLC right state
	std::array<double,4> UR_HLLC;
	UR_HLLC[0] = rho_R*((SR - vy_R)/(SR - Sstar));
	UR_HLLC[1] = UR_HLLC[0]*vx_R;
	UR_HLLC[2] = UR_HLLC[0]*Sstar;
	UR_HLLC[3] = UR_HLLC[0]*(E_R/rho_R + (Sstar - vy_R)*(Sstar + p_R/(rho_R*(SR - vy_R))));
	
	//Compute HLL flux 
	std::array<double, 4> fL_HLLC_y;
	std::array<double, 4> fR_HLLC_y; 
	for (int n = 0; n < 4; n++){
		fL_HLLC_y[n] = fy_L[n] + SL*(UL_HLLC[n] - UR_nph_y(i,j-1,k,n));
		fR_HLLC_y[n] = fy_R[n] + SR*(UR_HLLC[n] - UL_nph_y(i,j,k,n));
	}
	
	std::array<double, 4> finalflux;
	if (0 <= SL){
		finalflux = fy_L;
	}
	else if(SL < 0 && 0 <= Sstar){
		finalflux = fL_HLLC_y;
	}
	else if(Sstar < 0 &&  0 <= SR){
		finalflux = fR_HLLC_y;
	}
	else if(SR < 0){
		finalflux = fy_R; 
	}
		
	return finalflux;
}

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
  // Flux registers store fluxes at patch boundaries to ensure fluxes are conservative between AMR levels
  flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
  :
  AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
  // Flux registers are only required if AMR is actually used, and if flux fix up is being done (recommended)
  flux_reg = 0;
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
// AMReX can save simultion state such
// that if the code crashes, it can be restarted, with different
// settings files parameters if necessary (e.g. to output about the
// point of the crash).
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
  AmrLevel::restart(papa,is,bReadSpecial);
  
  BL_ASSERT(flux_reg == 0);
  if (level > 0 && do_reflux)
  {
    flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
  }
  
}

//
// Write a checkpoint file - format is handled automatically by AMReX
void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory - format is handled automatically by AMReX.
//
void
AmrLevelAdv::writePlotFile (const std::string& dir,
	 	            std::ostream&      os,
                            VisMF::How         how)
{
  AmrLevel::writePlotFile (dir,os,how);
}

//
//Define data descriptors.
//
// This is how the variables in a simulation are defined.  In the case
// of the advection equation, a single variable, phi, is defined.
//
void
AmrLevelAdv::variableSetUp ()
{
  BL_ASSERT(desc_lst.size() == 0);

  // A function which contains all processing of the settings file,
  // setting up initial data, choice of numerical methods and
  // boundary conditions
  read_params();
  
  const int storedGhostZones = 0;
    
  // Setting up a container for a variable, or vector of variables:
  // Phi_Type: Enumerator for this variable type
  // IndexType::TheCellType(): AMReX can support cell-centred and vertex-centred variables (cell centred here)
  // StateDescriptor::Point: Data can be a point in time, or an interval over time (point here)
  // storedGhostZones: Ghost zones can be stored (e.g. for output).  Generally set to zero.
  // NUM_STATE: Number of variables in the variable vector (1 in the case of advection equation)
  // cell_cons_interp: Controls interpolation between levels - cons_interp is good for finite volume
  desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
			 StateDescriptor::Point,storedGhostZones,NUM_STATE,
			 &cell_cons_interp);

  //Set up boundary conditions, all boundaries can be set
  //independently, including for individual variables, but lo (left) and hi (right) are useful ways to
  //store them, for consistent access notation for the boundary
  //locations
  int lo_bc[amrex::SpaceDim];
  int hi_bc[amrex::SpaceDim];
  // AMReX has pre-set BCs, including periodic (int_dir) and transmissive (foextrap)
  for (int i = 0; i < amrex::SpaceDim; ++i) {
    lo_bc[i] = hi_bc[i] = BCType::foextrap;   // periodic boundaries //need to change to transmissive for euler eqqations
  }

  // Object for storing all the boundary conditions
  BCRec bc(lo_bc, hi_bc);

  // Set up variable-specific information; needs to be done for each variable in NUM_STATE
  // Phi_Type: Enumerator for the variable type being set
  // 0: Position of the variable in the variable vector.  Single variable for advection.
  // phi: Name of the variable - appears in output to identify what is being plotted
  // bc: Boundary condition object for this variable (defined above)
  // BndryFunc: Function for setting boundary conditions.  For basic BCs, AMReX can handle these automatically
	
	// Added more components for the state vector - mc2246
	desc_lst.setComponent(Phi_Type, 0, "rho", bc, 
			StateDescriptor::BndryFunc(nullfill));
	desc_lst.setComponent(Phi_Type, 1, "momx", bc, 
			StateDescriptor::BndryFunc(nullfill));
	desc_lst.setComponent(Phi_Type, 2, "momy", bc, 
			StateDescriptor::BndryFunc(nullfill));
	desc_lst.setComponent(Phi_Type, 3, "Energy", bc, 
			StateDescriptor::BndryFunc(nullfill));						
	
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
  //
  // Loop over grids, call FORTRAN function to init with data.
  //
  const Real* dx  = geom.CellSize();
  // Position of the bottom left corner of the domain
  const Real* prob_lo = geom.ProbLo();
  // Create a multifab which can store the initial data
  MultiFab& S_new = get_new_data(Phi_Type);
  Real cur_time   = state[Phi_Type].curTime();

  // amrex::Print works like std::cout, but in parallel only prints from the root processor
  if (verbose) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  // Slightly messy way to ensure uninitialised data is not used.
  // AMReX has an XDim3 object, but a function needs to be written to
  // convert Real* to XDim3
  const Real dX = dx[0];
  // Is dimension greater than 1, if yes then use dy, if no then set to 0
  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
  // Is dimension greater than 2, if yes then use dz, if no then set to 0
  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);
  // get x0 domainMin
  const Real probLoX = prob_lo[0];
  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0); //get y0
  const Real probLoZ = (amrex::SpaceDim > 2 ? prob_lo[2] : 0.0); //get z0
  
  const Real gamma = 1.4;
  
  // Loop over all the patches at this level
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++){
      const Real z = probLoZ + (double(k)+0.5) * dZ;
      
      for(int j = lo.y; j <= hi.y; j++){
		const Real y = probLoY + (double(j)+0.5) * dY;
	
		for(int i = lo.x; i <= hi.x; i++){
			const Real x = probLoX + (double(i)+0.5) * dX;
	  
			if(amrex::SpaceDim == 1){ //Initial conditions in 1D
				const Real r2 =(x*x) / 0.01;
				arr(i,j,k,0) = 0.0 + exp(-r2);
				arr(i,j,k,1) = 1.0 + exp(-r2);
				arr(i,j,k,2) = 2.0 + exp(-r2);
			}
		
			else if(amrex::SpaceDim == 2){ //Initial conditions in 2D
				
				
				
				// Added code for setting up shock tube test - mc2246
				/*
				// Shock tube tests
				const Real rho_L = 1.0;
				const Real vx_L = 0.0;
				const Real vy_L = 0.0;
				const Real p_L = 1.0;
				
				const Real rho_R = 0.125;
				const Real vx_R = 0.0;
				const Real vy_R = 0.0;
				const Real p_R = 0.1;
				
				const Real x0 = 0.5;
				
				if (x < y){
				//if (y < x0){
				//if (x < x0){
					// Conservative variables
					arr(i,j,k,0) = rho_L; //rho
					arr(i,j,k,1) = rho_L*vx_L; //momx
					arr(i,j,k,2) = rho_L*vy_L; //momy
					arr(i,j,k,3) = p_L/(gamma - 1.0) + 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L); //E =p/(gamma -1) + 0.5*rho*(vx*vx + vy*vy)
				}
				
				else{
					// Conservative variables
					arr(i,j,k,0) = rho_R; //rho
					arr(i,j,k,1) = rho_R*vx_R; //momx
					arr(i,j,k,2) = rho_R*vy_R; //momy
					arr(i,j,k,3) = p_R/(gamma - 1.0) + 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R); //E =p/(gamma -1) + 0.5*rho*(vx*vx + vy*vy)
				}
				*/
				
				// Added code for setting up Toro Cylindrical Explosion test - mc2246
				// Toro Cylindrical Explosion
				const Real rho_L = 1.0;
				const Real vx_L = 0.0;
				const Real vy_L = 0.0;
				const Real p_L = 1.0;
				
				const Real rho_R = 0.125;
				const Real vx_R = 0.0; 
				const Real vy_R = 0.0;
				const Real p_R = 0.1;
				
				const Real R = 0.4;
				
				if ( ( (x-1.0)*(x-1.0) + (y-1.0)*(y-1.0) ) < R*R ){
					// Conservative variables
					arr(i,j,k,0) = rho_L; //rho
					arr(i,j,k,1) = rho_L*vx_L; //momx
					arr(i,j,k,2) = rho_L*vy_L; //momy
					arr(i,j,k,3) = p_L/(gamma - 1.0) + 0.5*rho_L*(vx_L*vx_L + vy_L*vy_L); //E =p/(gamma -1) + 0.5*rho*(vx*vx + vy*vy)
				}
				
				else{
					// Conservative variables
					arr(i,j,k,0) = rho_R; //rho
					arr(i,j,k,1) = rho_R*vx_R; //momx
					arr(i,j,k,2) = rho_R*vy_R; //momy
					arr(i,j,k,3) = p_R/(gamma - 1.0) + 0.5*rho_R*(vx_R*vx_R + vy_R*vy_R); //E =p/(gamma -1) + 0.5*rho*(vx*vx + vy*vy)
				}
				
				
			
			}
		
		  else{ //Initial conditions in 3D
			const Real r2 = (x*x + y*y + z*z) / 0.01;
			arr(i,j,k) = 1. + exp(-r2);
		  } 
		 
		  
		 }
		}
      }
    }

  if (verbose) {
    amrex::Print() << "Done initializing the level " << level 
		   << " data " << std::endl;
  }
}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init (AmrLevel &old)
{
  
  AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
  //
  // Create new grid data by fillpatching from old.
  //
  Real dt_new    = parent->dtLevel(level);
  Real cur_time  = oldlev->state[Phi_Type].curTime();
  Real prev_time = oldlev->state[Phi_Type].prevTime();
  Real dt_old    = cur_time - prev_time;
  setTimeLevel(cur_time,dt_old,dt_new);
  
  MultiFab& S_new = get_new_data(Phi_Type);

  const int zeroGhosts = 0;
  // FillPatch takes the data from the first argument (which contains
  // all patches at a refinement level) and fills (copies) the
  // appropriate data onto the patch specified by the second argument:
  // old: Source data
  // S_new: destination data
  // zeroGhosts: If this is non-zero, ghost zones could be filled too - not needed for init routines
  // cur_time: AMReX can attempt interpolation if a different time is specified - not recommended for advection eq.
  // Phi_Type: Specify the type of data being set
  // 0: This is the first data index that is to be copied
  // NUM_STATE: This is the number of states to be copied
  FillPatch(old, S_new, zeroGhosts, cur_time, Phi_Type, 0, NUM_STATE);

  // Note: In this example above, the all states in Phi_Type (which is
  // only 1 to start with) are being copied.  However, the FillPatch
  // command could be used to create a velocity vector from a
  // primitive variable vector.  In this case, the `0' argument is
  // replaced with the position of the first velocity component in the
  // primitive variable vector, and the NUM_STATE arguement with the
  // dimensionality - this argument is the number of variables that
  // are being filled/copied, and NOT the position of the final
  // component in e.g. the primitive variable vector.
}

//
//Initialize data on this level after regridding if old level did not previously exist
// These are standard AMReX commands which are unlikely to need altering
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);

    // See first init function for documentation
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//  This function is the one that actually calls the flux functions.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{

  MultiFab& S_mm = get_new_data(Phi_Type);

  // Note that some useful commands exist - the maximum and minumum
  // values on the current level can be computed directly - here the
  // max and min of variable 0 are being calculated, and output.
  Real maxval = S_mm.max(0);
  Real minval = S_mm.min(0);
  amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;

  // This ensures that all data computed last time step is moved from
  // `new' data to `old data' - this should not need changing If more
  // than one type of data were declared in variableSetUp(), then the
  // loop ensures that all of it is updated appropriately
  for (int k = 0; k < NUM_STATE_TYPE; k++) {
    state[k].allocOldData();
    state[k].swapTimeLevels(dt);
  }

  // S_new is the MultiFab that will be operated upon to update the data
  MultiFab& S_new = get_new_data(Phi_Type);

  const Real prev_time = state[Phi_Type].prevTime();
  const Real cur_time = state[Phi_Type].curTime();
  const Real ctr_time = 0.5*(prev_time + cur_time);

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

  //
  // Get pointers to Flux registers, or set pointer to zero if not there.
  //
  FluxRegister *fine    = 0; //need to set these to something else for the euler equations
  FluxRegister *current = 0;
    
  int finest_level = parent->finestLevel();

  // If we are not on the finest level, fluxes may need correcting
  // from those from finer levels.  To start this process, we set the
  // flux register values to zero
  if (do_reflux && level < finest_level) {
    fine = &getFluxReg(level+1);
    fine->setVal(0.0);
  }

  // If we are not on the coarsest level, the fluxes are going to be
  // used to correct those on coarser levels.  We get the appropriate
  // flux level to include our fluxes within
  if (do_reflux && level > 0)
  {
    current = &getFluxReg(level);
  }

  // Set up a dimensional multifab that will contain the fluxes
  MultiFab fluxes[amrex::SpaceDim];

  // Define the appropriate size for the flux MultiFab.
  // Fluxes are defined at cell faces - this is taken care of by the
  // surroundingNodes(j) command, ensuring the size of the flux
  // storage is increased by 1 cell in the direction of the flux.
  // This is only needed if refluxing is happening, otherwise fluxes
  // don't need to be stored, just used
  for (int j = 0; j < amrex::SpaceDim; j++)
  {
    BoxArray ba = S_new.boxArray();
    ba.surroundingNodes(j);
    fluxes[j].define(ba, dmap, NUM_STATE, 0);
  }

  // Advection velocity - AMReX allows the defintion of a vector
  // object (similar functionality to C++ std::array<N>, since its size must
  // be known, but was implemented before array was added to C++)
  //const Vector<Real> vel{1.0,1.0,0.0};
  
  // State with ghost cells - this is used to compute fluxes and perform the update.
  MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW);
  // See init function for details about the FillPatch function
  FillPatch(*this, Sborder, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
  
  Vector<BCRec> bcArray(NUM_STATE);
  
	// Modified to now use transmissive boundary conditions - mc2246
	for (int n = 0; n < NUM_STATE; ++n)
	{
		for (int idim = 0; idim < amrex::SpaceDim; ++idim)
		{
				bcArray[n].setLo(idim, BCType::foextrap); // first-order extrapolation
				bcArray[n].setHi(idim, BCType::foextrap);
		}
	}	
  
  // Fill periodic boundaries where they exist.  More accurately, the
  // FillBoundary call will fill overlapping boundaries (with periodic
  // domains effectively being overlapping).  It also takes care of
  // AMR patch and CPU boundaries.
  Sborder.FillBoundary(geom.periodicity());
  FillDomainBoundary(Sborder, geom, bcArray);
 
  
  

  for (int d = 0; d < amrex::SpaceDim ; d++)   
  {

    const int iOffset = ( d == 0 ? 1 : 0);
    const int jOffset = ( d == 1 ? 1 : 0);
    const int kOffset = ( d == 2 ? 1 : 0);
	
	
	if (d == 0){
		// Add code for x direction flux here
		
		
		// Added pragma loops of OpenMP - mc2246
		//****START PARALELL REGION HERE*****
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
		// Loop over all the patches at this level
		for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
		{
		  const Box& bx = mfi.tilebox();

		  const Dim3 lo = lbound(bx);
		  const Dim3 hi = ubound(bx);
		  
		  const Real* dx = geom.CellSize();
		  
		  const Real dX = dx[0];
		  // Is dimension greater than 1, if yes then use dy, if no then set to 0
		  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
		  // Is dimension greater than 2, if yes then use dz, if no then set to 0
		  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);
		  // get x0 domainMin

		  // Indexable arrays for the data, and the directional flux
		  // Based on the vertex-centred definition of the flux array, the
		  // data array runs from e.g. [0,N] and the flux array from [0,N+1]
		  
		  // Set up arrays to use in flux calculations - mc2246
		  const auto& arr = Sborder.array(mfi); //Conservative variables
		  const auto& UL_x = Sborder.array(mfi); //Left boundary extrapolated values
		  const auto& UR_x = Sborder.array(mfi); //Right boundary extrapolated values  
		  const auto& UL_nph_x = Sborder.array(mfi);
		  const auto& UR_nph_x = Sborder.array(mfi);
		   
		  const auto& fluxArr = fluxes[d].array(mfi);
		  		  
		  const Real omega  = 0.0; //put into inputs settings file
		  
		  
		  // Added loops for data reconstruction steps - mc2246
		  
		  // Need offset loops here because of indexing in x,y direction changes
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					for (int n = 0; n<4; n++){
						
						Real deltax = compute_deltax(arr,i,j,k,n,omega);
						Real r = compute_rx(arr,i,j,k,n);
						Real limiter = Minbeelimiter(r);
						
						UL_x(i,j,k,n) = arr(i,j,k,n) - 0.5*limiter*deltax;
						UR_x(i,j,k,n) = arr(i,j,k,n) + 0.5*limiter*deltax;
					}
				}
			}
		  }
		  
		  // Added loops for half time step evolution - mc2246
		  // Need offset loops here because of indexing in x,y direction changes
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					std::array<double,4> fluxL_x = feuler_x(UL_x,i,j,k);
					std::array<double,4> fluxR_x = feuler_x(UR_x,i,j,k);
					 
					for (int n = 0; n<4; n++){

						UL_nph_x(i,j,k,n) = UL_x(i,j,k,n) - 0.5*(dt/dX)*(fluxR_x[n] - fluxL_x[n]);
						UR_nph_x(i,j,k,n) = UR_x(i,j,k,n) - 0.5*(dt/dX)*(fluxR_x[n] - fluxL_x[n]);
						
						//UR_nph_x(i,j,k,n) = UR_x(i,j,k,n) + 0.5*(dt/dX)*(fluxR_x[n] - fluxL_x[n]);
					}
				}
			}
		  }
		  
		// Added loops for computing flux at each cell boundary interface - mc2246
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					//std::array<double,4> flux_x = feuler_x(arr,i,j,k);
					//std::array<double,4> flux_x = feuler_x(arr,i-iOffset,j-jOffset,k-kOffset);
					//std::array<double,4> flux_x = getHLLflux_x(UR_nph_x, UL_nph_x, i, j, k);
					std::array<double,4> flux_x = getHLLCflux_x(UR_nph_x, UL_nph_x, i, j, k);
					 
					// Insert HLLC flux code here 
					fluxArr(i,j,k,0) = flux_x[0];
					fluxArr(i,j,k,1) = flux_x[1];
					fluxArr(i,j,k,2) = flux_x[2];
					fluxArr(i,j,k,3) = flux_x[3];	
		
				}
			}
		  }

		  for(int k = lo.z; k <= hi.z; k++)
		  {
			for(int j = lo.y; j <= hi.y; j++)
			{
				for(int i = lo.x; i <= hi.x; i++)
				{
					for(int n = 0; n < 4; n++)
					{
						// Conservative update formula
						arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset,n) - fluxArr(i,j,k,n));
					
					}
				}
			}
		  }
		  
		}
	}
		//***END PARALLELL REGION HERE****

		// Modified to apply transmissive boundary conditions - mc2246
		// We need to compute boundary conditions again after each update
		Sborder.FillBoundary(geom.periodicity());
		FillDomainBoundary(Sborder, geom, bcArray);
		
		
		// The fluxes now need scaling for the reflux command.
		// This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
		if(do_reflux)
		{
		  Real scaleFactor = dt;
		  for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
		  {
		// Fluxes don't need scaling by dx[d]
		if(scaledir == d)
		{
		  continue;
		}
		scaleFactor *= dx[scaledir];
		  }
		  // The mult function automatically multiplies entries in a multifab by a scalar
		  // scaleFactor: The scalar to multiply by
		  // 0: The first data index in the multifab to multiply
		  // NUM_STATE:  The total number of data indices that will be multiplied
		  fluxes[d].mult(scaleFactor, 0, NUM_STATE);
		}
	 
	 
	 
	 
	} // close d==0 if statement
	
	
	
	else if (d == 1){
	 // Add code for y direction flux here
	 
	 // Added pragma loops of OpenMP - mc2246
	 //****START PARALELL REGION HERE*****
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
	 // Loop over all the patches at this level
		for (MFIter mfi(Sborder, true); mfi.isValid(); ++mfi)
		{
		  const Box& bx = mfi.tilebox();

		  const Dim3 lo = lbound(bx);
		  const Dim3 hi = ubound(bx);
		  
		  const Real* dx = geom.CellSize();
		  
		  const Real dX = dx[0];
		  // Is dimension greater than 1, if yes then use dy, if no then set to 0
		  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
		  // Is dimension greater than 2, if yes then use dz, if no then set to 0
		  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);
		  // get x0 domainMin

		  // Indexable arrays for the data, and the directional flux
		  // Based on the vertex-centred definition of the flux array, the
		  // data array runs from e.g. [0,N] and the flux array from [0,N+1]
		  
		  // Set up arrays to use in flux calculations - mc2246
		  const auto& arr = Sborder.array(mfi);
		  const auto& UL_y = Sborder.array(mfi); //Left boundary extrapolated values
		  const auto& UR_y = Sborder.array(mfi); //Right boundary extrapolated values  
		  const auto& UL_nph_y = Sborder.array(mfi); 
		  const auto& UR_nph_y = Sborder.array(mfi);
		  
		  const auto& fluxArr = fluxes[d].array(mfi);
		  
		  const Real omega  = 0.0; //put into inputs settings file
		  
		  
		  // Added loops for data reconstruction steps - mc2246
		   // Need offset loops here because of indexing in x,y direction changes
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					for (int n = 0; n<4; n++){
						
						Real deltay = compute_deltay(arr,i,j,k,n,omega);
						Real r = compute_ry(arr,i,j,k,n);
						Real limiter = Minbeelimiter(r);
						
						UL_y(i,j,k,n) = arr(i,j,k,n) - 0.5*limiter*deltay;
						UR_y(i,j,k,n) = arr(i,j,k,n) + 0.5*limiter*deltay;
					}
				}
			}
		  }
		  
		// Added loops for half time step evolution - mc2246
		  // Need offset loops here because of indexing in x,y direction changes
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					std::array<double,4> fluxL_y = feuler_y(UL_y,i,j,k);
					std::array<double,4> fluxR_y = feuler_y(UR_y,i,j,k);
					 
					for (int n = 0; n<4; n++){

						UL_nph_y(i,j,k,n) = UL_y(i,j,k,n) - 0.5*(dt/dY)*(fluxR_y[n] - fluxL_y[n]);
						UR_nph_y(i,j,k,n) = UR_y(i,j,k,n) - 0.5*(dt/dY)*(fluxR_y[n] - fluxL_y[n]);
						//UR_nph_y(i,j,k,n) = UR_y(i,j,k,n) + 0.5*(dt/dY)*(fluxR_y[n] - fluxL_y[n]);
					}
				}
			}
		  }
		  
		  // Added loops for computing flux at each cell boundary interface - mc2246
		  for(int k = lo.z; k <= hi.z+kOffset; k++)
		  {
			for(int j = lo.y; j <= hi.y+jOffset; j++)
			{
				for(int i = lo.x; i <= hi.x+iOffset; i++)
				{
					//std::array<double,4> flux_y = feuler_y(arr,i,j,k);
					//std::array<double,4> flux_y = feuler_y(arr,i-iOffset,j-jOffset,k-kOffset);
					//std::array<double,4> flux_y = getHLLflux_y(UR_nph_y, UL_nph_y, i, j, k);
					std::array<double,4> flux_y = getHLLCflux_y(UR_nph_y, UL_nph_y, i, j, k); 
					
					fluxArr(i,j,k,0) = flux_y[0];
					fluxArr(i,j,k,1) = flux_y[1];
					fluxArr(i,j,k,2) = flux_y[2];
					fluxArr(i,j,k,3) = flux_y[3];	
				
				}
			}
		  }
	
		  for(int k = lo.z; k <= hi.z; k++)
		  {
			for(int j = lo.y; j <= hi.y; j++)
			{
				for(int i = lo.x; i <= hi.x; i++)
				{
					for(int n = 0; n < 4; n++)
					{
					
						// Conservative update formula
						arr(i,j,k,n) = arr(i,j,k,n) - (dt / dx[d]) * (fluxArr(i+iOffset, j+jOffset, k+kOffset,n) - fluxArr(i,j,k,n));
					
					}
					
				}
			}
		  }
		  
		  
		}
	}
		//****END PARALELL REGION HERE****
		
		// Modified to apply transmissive boundary conditions - mc2246
		// We need to compute boundary conditions again after each update
		Sborder.FillBoundary(geom.periodicity());
		FillDomainBoundary(Sborder, geom, bcArray);
		
		
		// The fluxes now need scaling for the reflux command.
		// This scaling is by the size of the boundary through which the flux passes, e.g. the x-flux needs scaling by the dy, dz and dt
		if(do_reflux)
		{
		  Real scaleFactor = dt;
		  for(int scaledir = 0; scaledir < amrex::SpaceDim; ++scaledir)
		  {
		// Fluxes don't need scaling by dx[d]
		if(scaledir == d)
		{
		  continue;
		}
		scaleFactor *= dx[scaledir];
		  }
		  // The mult function automatically multiplies entries in a multifab by a scalar
		  // scaleFactor: The scalar to multiply by
		  // 0: The first data index in the multifab to multiply
		  // NUM_STATE:  The total number of data indices that will be multiplied
		  fluxes[d].mult(scaleFactor, 0, NUM_STATE);
		}
	 
	}
	
  }
  
  // The updated data is now copied to the S_new multifab.  This means
  // it is now accessible through the get_new_data command, and AMReX
  // can automatically interpolate or extrapolate between layers etc.
  // S_new: Destination
  // Sborder: Source
  // Third entry: Starting variable in the source array to be copied (the zeroth variable in this case)
  // Fourth entry: Starting variable in the destination array to receive the copy (again zeroth here)
  // NUM_STATE: Total number of variables being copied
  // Sixth entry: Number of ghost cells to be included in the copy (zero in this case, since only real
  //              data is needed for S_new)
  MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);

  // Refluxing at patch boundaries.  Amrex automatically does this
  // where needed, but you need to state a few things to make sure it
  // happens correctly:
  // FineAdd: If we are not on the coarsest level, the fluxes at this level will form part of the correction
  //          to a coarse level
  // CrseInit:  If we are not the finest level, the fluxes at patch boundaries need correcting.  Since we
  //            know that the coarse level happens first, we initialise the boundary fluxes through this
  //            function, and subsequently FineAdd will modify things ready for the correction
  // Both functions have the same arguments:
  // First: Name of the flux MultiFab (this is done dimension-by-dimension
  // Second: Direction, to ensure the correct vertices are being corrected
  // Third: Source component - the first entry of the flux MultiFab that is to be copied (it is possible that
  //        some variables will not need refluxing, or will be computed elsewhere (not in this example though)
  // Fourth: Destinatinon component - the first entry of the flux register that this call to FineAdd sends to
  // Fifth: NUM_STATE - number of states being added to the flux register
  // Sixth: Multiplier - in general, the least accurate (coarsest) flux is subtracted (-1) and the most
  //        accurate (finest) flux is added (+1)
  if (do_reflux) {
    if (current) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
    }
    if (fine) {
      for (int i = 0; i < amrex::SpaceDim ; i++)
	fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
    }
  }
  

  
  return dt;
}

//
//Estimate time step.
// This function is called by all of the other time step functions in AMReX, and is the only one that should
// need modifying
//
Real
AmrLevelAdv::estTimeStep (Real)
{
  // This is just a dummy value to start with 
  Real dt_est  = 1.0e+20;

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_new = get_new_data(Phi_Type);

  // This should not really be hard coded
  //const Real velMag = sqrt(2.);
  Real amax = 0.0;
  const Real gamma = 1.4;
  
  
  // Added pragmas for OpenMP - mc2246
  //**** START PARALELL REGION*****
  //need to do a reduction here 
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
  
  // Loop over all the patches at this level
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
  {
    Box bx = mfi.tilebox();
    const Dim3 lo = lbound(bx);
    const Dim3 hi = ubound(bx);

    const auto& arr = S_new.array(mfi);

    for(int k = lo.z; k <= hi.z; k++){
      //const Real z = probLoZ + (double(k)+0.5) * dZ;
      
      for(int j = lo.y; j <= hi.y; j++){
		//const Real y = probLoY + (double(j)+0.5) * dY;
	
		for(int i = lo.x; i <= hi.x; i++){
			//const Real x = probLoX + (double(i)+0.5) * dX;
			
			// Modified to compute timestep for Euler equations - mc2246
			
			Real rho = arr(i,j,k,0);
			Real momx = arr(i,j,k,1);
			Real momy = arr(i,j,k,2);
			Real E = arr(i,j,k,3);
			
			Real vx = arr(i,j,k,1)/arr(i,j,k,0);
			Real vy = arr(i,j,k,2)/arr(i,j,k,0);
			
			Real p = (gamma - 1.0)*(E - 0.5*rho*(vx*vx + vy*vy));
			
			Real Cs = sqrt(gamma*(p/rho));
			Real a = sqrt(vx*vx + vy*vy) + Cs;
			
			amax = std::max(amax, a); 
		 }
		
		}
		
      }
    
    }
     
  
  for(unsigned int d = 0; d < amrex::SpaceDim; ++d)
  {
    dt_est = std::min(dt_est, dx[d]/amax);
  }

  //****END PARALELL REGION HERE****
}
  
  // Ensure that we really do have the minimum across all processors
  ParallelDescriptor::ReduceRealMin(dt_est);
  dt_est *= cfl;

  if (verbose) {
    amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
		   << ":  dt_est = " << dt_est << std::endl;
  }
  
  return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
  return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
	  	               int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
  //
  // Grids have been constructed, compute dt for all levels.
  //
  // AMReX's AMR Level mode assumes that the time step only needs
  // calculating on the coarsest level - all subsequent time steps are
  // reduced by the refinement factor
  if (level > 0)
    return;

  // Initial guess
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor   *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_level[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
  //
  // We are at the end of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  //
  if (level > 0)
    return;

  // Although we only compute the time step on the finest level, we
  // need to take information from all levels into account.  The
  // sharpest features may be smeared out on coarse levels, so not
  // using finer levels could cause instability
  for (int i = 0; i <= finest_level; i++)
  {
    AmrLevelAdv& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  // A couple of things are implemented to ensure that time step's
  // don't suddenly grow by a lot, as this could lead to errors - for
  // sensible mesh refinement choices, these shouldn't really change
  // anything
  if (post_regrid_flag == 1) 
  {
    //
    // Limit dt's by pre-regrid dt
    //
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],dt_level[i]);
    }
  }
  else 
  {
    //
    // Limit dt's by change_max * old dt
    //
    static Real change_max = 1.1;
    for (int i = 0; i <= finest_level; i++)
    {
      dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
    }
  }
    
  //
  // Find the minimum over all levels
  //
  Real dt_0 = 1.0e+100;
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_0 = std::min(dt_0,n_factor*dt_min[i]);
  }

  //
  // Limit dt's by the value of stop_time.
  //
  const Real eps = 0.001*dt_0;
  Real cur_time  = state[Phi_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - eps))
      dt_0 = stop_time - cur_time;
  }
  
  n_factor = 1;
  for (int i = 0; i <= finest_level; i++)
  {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0/n_factor;
  }
}

//
//Do work after timestep().
// If something has to wait until all processors have done their advance function, the post_timestep function
// is the place to put it.  Refluxing and averaging down ar            e two standard examples for AMR
//
void
AmrLevelAdv::post_timestep (int iteration)
{
  //
  // Integration cycle on fine level grids is complete
  // do post_timestep stuff here.
  //
  
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const Real cur_time = state[Phi_Type].curTime();
  const MultiFab& S_new = get_new_data(Phi_Type);
  
  
  
  
  int finest_level = parent->finestLevel();
  
  if (do_reflux && level < finest_level)
    reflux();
  
  if (level < finest_level)
    avgDown();
  // Added function which outputs data from each level  into a .txt file at the time 'finaltime' - mc2246
  //Move to post timestep function
  Real finaltime = 0.25;
  
  if (cur_time == finaltime){
	  
	  std::string levelnum = std::to_string(level);
	  
	  std::ofstream outputfinal("finalplots_"+ levelnum + ".txt");
	  
	  const Real* dx  = geom.CellSize();
	  const Real* prob_lo = geom.ProbLo();

	  const Real dX = dx[0];
	  const Real dY = (amrex::SpaceDim > 1 ? dx[1] : 0.0);
	  const Real dZ = (amrex::SpaceDim > 2 ? dx[2] : 0.0);
	  
	  const Real probLoX = prob_lo[0];
	  const Real probLoY = (amrex::SpaceDim > 1 ? prob_lo[1] : 0.0); //get y0
	  const Real probLoZ = (amrex::SpaceDim > 2 ? prob_lo[2] : 0.0); //get z0
	  
	  
	  const Real gamma = 1.4;
	  
	  // Loop over all the patches at this level
	  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	  {
		Box bx = mfi.tilebox();
		const Dim3 lo = lbound(bx);
		const Dim3 hi = ubound(bx);

		const auto& arr = S_new.array(mfi);

		for(int k = lo.z; k <= hi.z; k++){
		  const Real z = probLoZ + (double(k)+0.5) * dZ;
		  
		  for(int j = lo.y; j <= hi.y; j++){
				const Real y = probLoY + (double(j)+0.5) * dY;
			
			for(int i = lo.x; i <= hi.x; i++){
				const Real x = probLoX + (double(i)+0.5) * dX;

					const Real gamma = 1.4;
					
					const Real rho = arr(i,j,k,0);
					const Real momx = arr(i,j,k,1);
					const Real momy = arr(i,j,k,2);
					const Real E = arr(i,j,k,3);
					
					const Real vx = momx/rho;
					const Real vy = momy/rho;
					const Real p = (gamma - 1.0)*(E - 0.5*rho*(vx*vx + vy*vy));
					const Real internalenergy =  p/(rho*(gamma - 1.0));
					
					if (y == (probLoY + (double(lo.y)+0.5) * dY)){
						//outputfinal << x << " " << y << " " << rho << " " << vx << " " << vy << " " << p << " " << internalenergy << " " << level << std::endl;
						outputfinal << x << " " << rho << std::endl;
					}
			 }
			 //Add extra blank line
			//outputfinal << " " << std::endl;
			 
			}
		  }
	   }
    }
    
    
  
}

//
//Do work after regrid().
// Nothing normally needs doing here, but if something was calculated on a per-patch basis, new patches might
// this to be calcuated immediately
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest)
{

}

//
//Do work after a restart().
// Similar to post_regrid, nothing normally needs doing here
//
void
AmrLevelAdv::post_restart() 
{

}

//
//Do work after init().
// Once new patches have been initialised, work may need to be done to ensure consistency, for example,
// averaging down - though for linear interpolation, this probably won't change anything
//
void
AmrLevelAdv::post_init (Real stop_time)
{
  if (level > 0)
    return;
  //
  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  //
  int finest_level = parent->finestLevel();
  for (int k = finest_level-1; k>= 0; k--)
    getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//  Determine which parts of the domain need refinement
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
	
   amrex::Print() << "Doing error estimation on level " << level 
		   << " data " << std::endl;
		   
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  MultiFab& S_new = get_new_data(Phi_Type);

	//***START PARALELL REGION HERE****
//#ifdef AMREX_USE_OMP
//#pragma omp parallel
//#endif
    //{

  Vector<int> itags;
	
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
  {
    const Box&  tilebx  = mfi.tilebox();

    // An AMReX construction, effectively a boolean array which is true in positions that are valid for refinement
    TagBox&     tagfab  = tags[mfi];

    // Traditionally, a lot of the array-based operations in AMReX happened in Fortran.  The standard template
    // for these is short and easy to read, flagging on values or gradients (first order calculation)
    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
    // So we are going to get a temporary integer array.
    tagfab.get_itags(itags, tilebx);
	    
    // data pointer and index space
    int*        tptr    = itags.dataPtr();
    const int*  tlo     = tilebx.loVect();
    const int*  thi     = tilebx.hiVect();

    // Various macros exist to convert the C++ data structures to Fortran
    state_error(tptr,  AMREX_ARLIM_3D(tlo), AMREX_ARLIM_3D(thi),
		BL_TO_FORTRAN_3D(S_new[mfi]),
		&tagval, &clearval, 
		AMREX_ARLIM_3D(tilebx.loVect()), AMREX_ARLIM_3D(tilebx.hiVect()), 
		AMREX_ZFILL(dx), AMREX_ZFILL(prob_lo), &time, &level);
    //
    // Now update the tags in the TagBox.
    //
    tagfab.tags_and_untags(itags, tilebx);
  }
  
  //***END PARALELL REGION HERE****
	//}
}

//
// This function reads the settings file
//
void
AmrLevelAdv::read_params ()
{
  // Make sure that this is only done once
  static bool done = false;

  if (done) return;

  done = true;

  // A ParmParse object allows settings, with the correct prefix, to be read in from the settings file
  // The prefix can help identify what a settings parameter is used for
  // AMReX has some default ParmParse names, amr and geometry are two commonly needed ones
  ParmParse pp("adv");   

  // ParmParse has two options; query and get.  Query will only alter
  // a parameter if it can be found (if these aren't in the settings
  // file, then the values at the top of this file will be used).  Get
  // will throw an error if the parameter is not found in the settings
  // file.
  pp.query("v",verbose);
  pp.query("cfl",cfl);
  pp.query("do_reflux",do_reflux);

  // Vector variables can be read in; these require e.g.\ pp.queryarr
  // and pp.getarr, so that the ParmParse object knows to look for
  // more than one variable

  // Geometries can be Cartesian, cylindrical or spherical - some
  // functions (e.g. divergence in linear solvers) are coded with this
  // geometric dependency
  Geometry const* gg = AMReX::top()->getDefaultGeometry();

  // This tutorial code only supports Cartesian coordinates.
  if (! gg->IsCartesian()) {
    amrex::Abort("Please set geom.coord_sys = 0");
  }

  // This tutorial code only supports periodic boundaries.
  // The periodicity is read from the settings file in AMReX source code, but can be accessed here
  //if (! gg->isAllPeriodic()) {
  //  amrex::Abort("Please set geometry.is_periodic = 1 1 1");
  //}

  //
  // read tagging parameters from probin file
  //
  // Tradtionally, the inputs file with ParmParse functionality is handled by C++.  However, a Fortran settings
  // file, by default named probin, can also supply variables.  Mostly used for mesh refinement (tagging) critera
  std::string probin_file("probin");

  ParmParse ppa("amr");
  ppa.query("probin_file",probin_file);

  int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  // use a fortran routine to
  // read in tagging parameters from probin file
  get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);

}

//
// AMReX has an inbuilt reflux command, but we still have the freedom
// to decide what goes into it (for example, which variables are
// actually refluxed).  This also gives a little flexibility as to
// where flux registers are stored.  In this example, they are stored
// on levels [1,fine] but not level 0.  
//
void
AmrLevelAdv::reflux ()
{
  BL_ASSERT(level<parent->finestLevel());

  const Real strt = amrex::second();

  // Call the reflux command with the appropriate data.  Because there
  // are no flux registers on the coarse level, they start from the
  // first level.  But the coarse level to the (n-1)^th are the ones
  // that need refluxing, hence the `level+1'.  
  getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
  if (verbose)
  {
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    Real      end    = amrex::second() - strt;
    
    ParallelDescriptor::ReduceRealMax(end,IOProc);
    
    amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
		   << " : time = " << end << std::endl;
  }
}

//
// Generic function for averaging down - in this case it just makes sure it doesn't happen on the finest level
//
void
AmrLevelAdv::avgDown ()
{
  if (level == parent->finestLevel())
  {
    return;
  }
  // Can select which variables averaging down will happen on - only one to choose from in this case!
  avgDown(Phi_Type);
}

//
// Setting up the call to the AMReX-implemented average down function
//
void
AmrLevelAdv::avgDown (int state_indx)
{
  // For safety, again make sure this only happens if a finer level exists
  if (level == parent->finestLevel()) return;

  // You can access data at other refinement levels, use this to
  // specify your current data, and the finer data that is to be
  // averaged down
  AmrLevelAdv& fine_lev = getLevel(level+1);
  MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
  MultiFab&  S_crse   = get_new_data(state_indx);

  // Call the AMReX average down function:
  // S_fine: Multifab with the fine data to be averaged down
  // S_crse: Multifab with the coarse data to receive the fine data where necessary
  // fine_lev.geom:  Geometric information (cell size etc.) for the fine level
  // geom: Geometric information for the coarse level (i.e. this level)
  // 0: First variable to be averaged (as not all variables need averaging down
  // S_fine.nComp(): Number of variables to average - this can be computed automatically from a multifab
  // refRatio: The refinement ratio between this level and the finer level
  amrex::average_down(S_fine,S_crse,
		      fine_lev.geom,geom,
		      0,S_fine.nComp(),parent->refRatio(level));
}
