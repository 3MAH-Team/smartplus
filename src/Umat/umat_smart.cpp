///@file abaqus2smart.hpp
///@brief Procedure that transfer the abaqus format to a SMART+ format:
///@brief Implemented in 1D-2D-3D
///@author Chemisky
///@version 1.0
///@date 12/03/2013

#include <iostream>
#include <map>
#include <fstream>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <armadillo>
#include <smartplus/parameter.hpp>
#include <smartplus/Umat/umat_smart.hpp>

#include <smartplus/Umat/Mechanical/Elasticity/elastic_isotropic.hpp>
#include <smartplus/Umat/Mechanical/Elasticity/elastic_transverse_isotropic.hpp>
#include <smartplus/Umat/Mechanical/Plasticity/plastic_isotropic_ccp.hpp>

#include <smartplus/Umat/Thermomechanical/Elasticity/elastic_isotropic.hpp>

#include <smartplus/Micromechanics/Mori_Tanaka/Mori_Tanaka.hpp>
#include <smartplus/Micromechanics/Periodic_Layer/Periodic_Layer.hpp>
#include <smartplus/Micromechanics/Self_Consistent/Self_Consistent.hpp>

using namespace std;
using namespace arma;

namespace smart{

///@param stress array containing the components of the stress tensor (dimension ntens)
///@param stran array containing total strain component (dimension ntens) at the beginning of increment
///@param dstran array containing the component of total strain increment (dimension ntens)
///@param time two compoenent array : first component is the value of step time at the beginning of the current increment and second component is the value of total time at the beginning of the current increment
///@param dtime time increment
///@param temperature temperature avlue at the beginning of increment
///@param Dtemperature temperature increment
///@param ndi number of direct stress components
///@param nshr number of shear stress components
///@param drot rotation increment matrix (dimension 3*3)
///@param pnewdt ratio of suggested new time increment
///@param celent characteristic element length
///@param dfgrd0 array containing the deformation gradient at the beginning of increment (dimension 3*3)
///@param dfgrd1 array containing the deformation gradient at the end of increment (dimension 3*3)
///@param noel element number
///@param npt integration point number
///@param layer layer number - not used
///@param kspt section point number within the current layer - not used
///@param kstep step number
///@param kinc increment number


void abaqus2smart(double *stress, double *ddsdde, const double *stran, const double *dstran, const double *time, const double &dtime, const double &temperature, const double &Dtemperature, const int &nprops,const double *props, const int &nstatev, double *statev, const double &pnewdt, const int &ndi, const int &nshr, const double *drot, vec &sigma, mat &Lt, vec &Etot, vec &DEtot, double &T, double &DT, double &Time, double &DTime, vec &props_smart, vec &statev_smart, double &tnew_dt, mat &DR, bool &start)
{
    
	if(ndi == 1){						// 1D
		sigma(0) = stress[0];
		Etot(0) = stran[0];
		DEtot(0) = dstran[0];
		Lt(0,0) = ddsdde[0];
		
	}
	else if(ndi == 2){					// 2D Plane Stress
		sigma(0) = stress[0];
		sigma(1) = stress[1];
		sigma(3) = stress[2];
		
		Etot(0) = stran[0];
		Etot(1) = stran[1];
		Etot(3) = stran[2];
        
		DEtot(0) = dstran[0];
		DEtot(1) = dstran[1];
		DEtot(3) = dstran[2];
		      
        for(int i=0 ; i<3 ; i++)
        {
            for(int j=0 ; j<3 ; j++)
				Lt(j,i) = ddsdde[i*3+j];
        }
	}
	else if(ndi == 3){
		if(nshr == 1) {
			sigma(0) = stress[0];		// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
			sigma(1) = stress[1];
			sigma(2) = stress[2];
			sigma(3) = stress[3];
			
			Etot(0) = stran[0];
			Etot(1) = stran[1];
			Etot(2) = stran[2];
			Etot(3) = stran[3];
            
			DEtot(0) = dstran[0];
			DEtot(1) = dstran[1];
			DEtot(2) = dstran[2];
			DEtot(3) = dstran[3];
			
			for(int i=0 ; i<4 ; i++)
			{
				for(int j=0 ; j<4 ; j++)
				Lt(j,i) = ddsdde[i*4+j];
			}							
			
		}
		else {							// 3D
			sigma(0) = stress[0];
			sigma(1) = stress[1];
			sigma(2) = stress[2];
			sigma(3) = stress[3];
			sigma(4) = stress[4];
			sigma(5) = stress[5];
			
			Etot(0) = stran[0];
			Etot(1) = stran[1];
			Etot(2) = stran[2];
			Etot(3) = stran[3];
			Etot(4) = stran[4];
			Etot(5) = stran[5];
            
			DEtot(0) = dstran[0];
			DEtot(1) = dstran[1];
			DEtot(2) = dstran[2];
			DEtot(3) = dstran[3];
			DEtot(4) = dstran[4];
			DEtot(5) = dstran[5];

			for(int i=0 ; i<6 ; i++)
			{
				for(int j=0 ; j<6 ; j++)
				Lt(j,i) = ddsdde[i*6+j];
			}				
			
		}
	}
	
	///@brief rotation matrix
	DR(0,0) = drot[0];
	DR(0,1) = drot[3];
	DR(0,2) = drot[6];
	DR(1,0) = drot[1];
	DR(1,1) = drot[4];
	DR(1,2) = drot[7];
	DR(2,0) = drot[2];
	DR(2,1) = drot[5];
	DR(2,2) = drot[8];
	
	///@brief Temperature and temperature increment creation
	T = temperature;
	DT = Dtemperature;
    tnew_dt = pnewdt;
    
	///@brief Time
	Time = time[1];
	DTime = dtime;
    
	///@brief Initialization
	if(Time < 1E-12)
	{
		start = true;
	}
	else
	{
		start = false;
	}
    
    ///@brief : Pass the material properties and the internal variables
    for (int i=0; i<nprops; i++) {
        props_smart(i) = props[i];
    }
    for (int i=0; i<nstatev; i++) {
        statev_smart(i) = statev[i];
    }
	
}

void select_umat_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &rpl, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, bool &start, double &tnew_dt)
{
    string buffer;
    int c_umat = 0;
    int n_umat;
    int n_meta;
    
    ifstream listumat;
    listumat.open("data/list_umat_T.dat", ios::in);
    if(listumat) {
        listumat >> buffer >> n_umat >> buffer >> n_meta;
    }
    else {
        cout << "Error: cannot open list_umat_T.dat file \n";
    }
    
    std::map<string, int> list_umat;
    for (int i = 0; i < (n_umat + n_meta); i++){
        listumat >> buffer >> buffer >> c_umat;
        list_umat.insert ( std::pair<string,int>(buffer,c_umat) );
    }
    listumat.close();
    
    switch (list_umat[umat_name]) {
        case 1: {
            umat_elasticity_iso_T(Etot, DEtot, sigma, rpl, dSdE, dSdT, drpldE, drpldT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
            break;
        }
            
        default: {
            cout << "Error: The choice of Thermomechanical Umat could not be found in the umat library :" << umat_name << "\n"; 
        }
    }
    
}
    
void select_umat(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, const bool &start, double &tnew_dt)
{
	string buffer;
    int c_umat = 0;
	int n_umat;
	int n_meta;
	
	ifstream listumat;
	listumat.open("data/list_umat.dat", ios::in);
	if(listumat) {
		listumat >> buffer >> n_umat >> buffer >> n_meta;
	}
	else {
		cout << "Error: cannot open list_umat.dat file \n";
	}
	
    std::map<string, int> list_umat;
    for (int i = 0; i < (n_umat + n_meta); i++){
		listumat >> buffer >> buffer >> c_umat;
        list_umat.insert ( std::pair<string,int>(buffer,c_umat) );
    }
	listumat.close();
    
	switch (list_umat[umat_name]) {
		case 1: {
			umat_elasticity_iso(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
		case 2: {
			umat_elasticity_trans_iso(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
		case 3: {
			umat_plasticity_iso_CCP(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
		case 100: {
			umat_MT_N(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
		case 101: {
			umat_SC_N(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
		case 102: {
			umat_PL_N(Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start);
			break;
		}
						
		default: {
			cout << "Error: The choice of Umat could not be found in the umat library :" << umat_name << "\n"; 
		}
	}
    
}

void run_umat_T(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, double &rpl, mat &dSdE, mat &dSdT, mat &drpldE, mat &drpldT, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, bool &start, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    select_umat_T(umat_name, Etot, DEtot, sigma, rpl, dSdE, dSdT, drpldE, drpldT, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
    
    if (Time + DTime > precision) {
        start = false;
    }
}

void run_umat(const string &umat_name, const vec &Etot, const vec &DEtot, vec &sigma, mat &Lt, const mat &DR, const int &nprops, const vec &props, const int &nstatev, vec &statev, const double &T, const double &DT,const double &Time,const double &DTime, double &sse, double &spd, const int &ndi, const int &nshr, bool &start, double &tnew_dt)
{
    
    tnew_dt = 1.;
    
    select_umat(umat_name, Etot, DEtot, sigma, Lt, DR, nprops, props, nstatev, statev, T, DT, Time, DTime, sse, spd, ndi, nshr, start, tnew_dt);
    
    if (Time + DTime > precision) {
        start = false;
    }
}	

void smart2abaqus(double *stress, double *ddsdde, double *statev, const int &ndi, const int &nshr, const vec &sigma, const mat &Lt, const vec &statev_smart, double &pnewdt, const double &tnew_dt)
{
 
    pnewdt = tnew_dt;
    
    if(ndi == 1) {							// 1D
        stress[0] = sigma(0);
        ddsdde[0] = Lt(0,0);
    }
    else if(ndi == 2) {						// 2D Plane Stress
        stress[0] = sigma(0);
        stress[1] = sigma(1);
        stress[2] = sigma(3);
        
        assert(Lt(2,2) > 0.);
        ddsdde[0] = Lt(0,0)-Lt(0,2)*Lt(2,0)/Lt(2,2);
        ddsdde[4] = Lt(1,1)-Lt(1,2)*Lt(2,1)/Lt(2,2);
        ddsdde[3] = Lt(0,1)-Lt(0,2)*Lt(2,1)/Lt(2,2);
        ddsdde[1] = Lt(1,0)-Lt(1,2)*Lt(2,0)/Lt(2,2);
        ddsdde[6] = Lt(0,3)-Lt(0,2)*Lt(2,3)/Lt(2,2);
        ddsdde[7] = Lt(1,3)-Lt(1,2)*Lt(2,3)/Lt(2,2);
        ddsdde[2] = Lt(3,0)-Lt(3,2)*Lt(2,0)/Lt(2,2);
        ddsdde[5] = Lt(3,1)-Lt(3,2)*Lt(2,1)/Lt(2,2);
        ddsdde[8] = Lt(3,3)-Lt(3,2)*Lt(2,3)/Lt(2,2);
    }
    else if(ndi == 3){
		if (nshr == 1) {					// 2D Generalized Plane Strain (Plane Strain, Axisymetric)
			stress[0] = sigma(0);
			stress[1] = sigma(1);
			stress[2] = sigma(2);
			stress[3] = sigma(3);
			
			ddsdde[0] = Lt(0,0);
			ddsdde[4] = Lt(0,1);
			ddsdde[8] = Lt(0,2);
			ddsdde[3] = Lt(0,3);
			ddsdde[1] = Lt(1,0);
			ddsdde[5] = Lt(1,1);
			ddsdde[9] = Lt(1,2);
			ddsdde[7] = Lt(1,3);
			ddsdde[2] = Lt(2,0);
			ddsdde[6] = Lt(2,1);
			ddsdde[10] = Lt(2,2);
			ddsdde[11] = Lt(2,3);
			ddsdde[12] = Lt(3,0);
			ddsdde[13] = Lt(3,1);
			ddsdde[14] = Lt(3,2);
			ddsdde[15] = Lt(3,3);
		}
		else {								// 3D
			stress[0] = sigma(0);
			stress[1] = sigma(1);
			stress[2] = sigma(2);
			stress[3] = sigma(3);
			stress[4] = sigma(4);
			stress[5] = sigma(5);
			
			ddsdde[0] = Lt(0,0);
			ddsdde[6] = Lt(0,1);
			ddsdde[12] = Lt(0,2);
			ddsdde[18] = Lt(0,3);
			ddsdde[24] = Lt(0,4);
			ddsdde[30] = Lt(0,5);
			ddsdde[1] = Lt(1,0);
			ddsdde[7] = Lt(1,1);
			ddsdde[13] = Lt(1,2);
			ddsdde[19] = Lt(1,3);
			ddsdde[25] = Lt(1,4);
			ddsdde[31] = Lt(1,5);
			ddsdde[2] = Lt(2,0);
			ddsdde[8] = Lt(2,1);
			ddsdde[14] = Lt(2,2);
			ddsdde[20] = Lt(2,3);
			ddsdde[26] = Lt(2,4);
			ddsdde[32] = Lt(2,5);
			ddsdde[3] = Lt(3,0);
			ddsdde[9] = Lt(3,1);
			ddsdde[15] = Lt(3,2);
			ddsdde[21] = Lt(3,3);
			ddsdde[27] = Lt(3,4);
			ddsdde[33] = Lt(3,5);
			ddsdde[4] = Lt(4,0);
			ddsdde[10] = Lt(4,1);
			ddsdde[16] = Lt(4,2);
			ddsdde[22] = Lt(4,3);
			ddsdde[28] = Lt(4,4);
			ddsdde[34] = Lt(4,5);
			ddsdde[5] = Lt(5,0);
			ddsdde[11] = Lt(5,1);
			ddsdde[17] = Lt(5,2);
			ddsdde[23] = Lt(5,3);
			ddsdde[29] = Lt(5,4);
			ddsdde[35] = Lt(5,5);
		}
	}
	
    ///@brief : Pass the material properties and the variables
    for (int i=0; i<statev_smart.n_elem; i++) {
        statev[i] = statev_smart(i);
    }
    
    
}
	
} //namespace smart
