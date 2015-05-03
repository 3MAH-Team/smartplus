/* This file is part of SMART+.
 
 SMART+ is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 SMART+ is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with SMART+.  If not, see <http://www.gnu.org/licenses/>.
 
 */

///@file ellipsoid_characteristics.cpp
///@brief characteristics of an ellipsoid
///@version 1.0

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <armadillo>
#include <smartplus/Libraries/Homogenization/phase_characteristics.hpp>
#include <smartplus/Libraries/Homogenization/ellipsoid_characteristics.hpp>
#include <smartplus/Libraries/Maths/rotation.hpp>
#include <smartplus/Libraries/Homogenization/eshelby.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================

/*!
  \brief default constructor
*/

//-------------------------------------------------------------
ellipsoid_characteristics::ellipsoid_characteristics() : phase_characteristics(), S(6,6), T(6,6)
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
    
	a1=0.;
	a2=0.;
    a3=0.;
	
    psi_geom=0.;
    theta_geom=0.;
    phi_geom=0.;
    
	S = zeros(6,6);
	T = zeros(6,6);
}

/*!
  \brief Constructor
  \param nprops : size of the table statev 
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
ellipsoid_characteristics::ellipsoid_characteristics(int n, int m, bool init, double value) : phase_characteristics(n,m,init,value), S(6,6), T(6,6)
//-------------------------------------------------------------
{
    coatingof = 0;
    coatedby = 0;
    
	a1=0.;
	a2=0.;
    a3=0.;
	
    psi_geom=0.;
    theta_geom=0.;
    phi_geom=0.;
    
	S = zeros(6,6);
	T = zeros(6,6);
}

/*!
  \brief Constructor with parameters
  \param nstatev : size of the table statev 
  \param init boolean that indicates if the constructor has to initialize the state variables vectors and statev table (default value is true) \n
  \n\n
  \f$ \textbf{Examples :} \f$ \n
*/

//-------------------------------------------------------------
ellipsoid_characteristics::ellipsoid_characteristics(int mnumber, int mcoatingof, int mcoatedby, string mumat_name, double mconcentration, double mpsi_mat, double mtheta_mat, double mphi_mat, double mpsi_geom, double mtheta_geom, double mphi_geom, double ma1, double ma2, double ma3, int mnprops, const vec &mprops, int mnstatev, const vec &mstatev, const state_variables &mlocal, const state_variables &mglobal, const mat &mS, const mat &mT, const mat &mA, const mat &mB) : phase_characteristics(mnumber, mumat_name, mconcentration, mpsi_mat, mtheta_mat, mphi_mat, mnprops, mprops, mnstatev, mstatev, mlocal, mglobal, mA, mB), S(6,6), T(6,6)
//-------------------------------------------------------------
{	

	assert (mS.n_rows == 6);
	assert (mS.n_cols == 6);
	assert (mT.n_rows == 6);
	assert (mT.n_cols == 6);
	
	coatingof = mcoatingof;
	coatedby = mcoatedby;
    
	psi_geom = mpsi_geom;
	theta_geom = mtheta_geom;
	phi_geom = mphi_geom;
	a1 = ma1;
	a2 = ma2;
	a3 = ma3;
	
	S = mS;
	T = mT;
}

/*!
  \brief Copy constructor
  \param s phases_characteristics object to duplicate
*/

//------------------------------------------------------
ellipsoid_characteristics::ellipsoid_characteristics(const ellipsoid_characteristics& ev) : phase_characteristics(ev), S(6,6), T(6,6)
//------------------------------------------------------
{
	coatingof = ev.coatingof;
	coatedby = ev.coatedby;
    
	psi_geom = ev.psi_geom;
	theta_geom = ev.theta_geom;
	phi_geom = ev.phi_geom;
	a1 = ev.a1;
	a2 = ev.a2;
	a3 = ev.a3;
	
	S = ev.S;
	T = ev.T;
}

/*!
  \brief Destructor

  Deletes statev_variables, the vectors and matrix, and table statev if nstatev is not null.
*/

//-------------------------------------
ellipsoid_characteristics::~ellipsoid_characteristics()
//-------------------------------------
{ }

void ellipsoid_characteristics::fillS(const mat& L_m, const vec &x, const vec &wx, const vec &y, const vec &wy, const int &mp, const int &np)
{
    mat Lm_local_geom = rotate_g2l_L(L_m, psi_geom, theta_geom, phi_geom);
    S = Eshelby(Lm_local_geom, a1, a2, a3, x, wx, y, wy, mp, np);
}

void ellipsoid_characteristics::fillT(const mat& L_m, const vec &x, const vec &wx, const vec &y, const vec &wy, const int &mp, const int &np)
{
    mat Lm_local_geom = rotate_g2l_L(L_m, psi_geom, theta_geom, phi_geom);
    S = Eshelby(Lm_local_geom, a1, a2, a3, x, wx, y, wy, mp, np);
    mat Lt_local_geom = rotate_g2l_L(global.Lt, psi_geom, theta_geom, phi_geom);
    
    mat T_local_geom = inv(eye(6,6) + S*inv(Lm_local_geom)*(Lt_local_geom - Lm_local_geom));
    
    T = rotate_l2g_A(T_local_geom, psi_geom, theta_geom, phi_geom);
}

/*!
 \brief Standard operator = for phases_characteristics
 */
//----------------------------------------------------------------------
ellipsoid_characteristics& ellipsoid_characteristics::operator = (const ellipsoid_characteristics& ev)
//----------------------------------------------------------------------
{
	assert(ev.nprops);

	number = ev.number;
	coatingof = ev.coatingof;
    coatedby = ev.coatedby;
	umat_name = ev.umat_name;

	concentration = ev.concentration;
	psi_mat = ev.psi_mat;
	theta_mat = ev.theta_mat;
	phi_mat = ev.phi_mat;
    
	psi_geom = ev.psi_geom;
	theta_geom = ev.theta_geom;
	phi_geom = ev.phi_geom;
    
	a1 = ev.a1;
	a2 = ev.a2;
	a3 = ev.a3;

	local = ev.local;
	global = ev.global;
	
	S = ev.S;
	T = ev.T;
	A = ev.A;
	B = ev.B;
	
	nprops = ev.nprops;
	props = ev.props;
	nstatev = ev.nstatev;
	statev = ev.statev;

	return *this;

}

//--------------------------------------------------------------------------
ostream& operator << (ostream& s, const ellipsoid_characteristics& ev)
//--------------------------------------------------------------------------
{
	assert(ev.nprops);

	s << "Display state variables\n";
	s << "Number of the phase: " << ev.number << "\n";
	s << "is a coating of the phase: " << ev.coatingof << "\n";
	s << "is coated by the phase: " << ev.coatedby << "\n";
	s << "Name of the material = " << ev.umat_name << "\n";
	s << "concentration: " << ev.concentration << "\n";
	s << "local material orientation: psi = " << ev.psi_mat << "\t theta = " << ev.theta_mat << "\t phi = " << ev.phi_mat << "\n";
	s << "local geometrical orientation: psi = " << ev.psi_geom << "\t theta = " << ev.theta_geom << "\t phi = " << ev.phi_geom << "\n";
	s << "geometry: " << ev.a1 << "\t" << ev.a2 << "\t" << ev.a3 << "\n";
	s << "nprops: \n" << ev.nprops << "\n";
	s << "props: \n";
    s << ev.props.t();
    s << "\n";

	s << "local state variables: \n" << ev.local << "\n";
	s << "global state variables: \n" << ev.global << "\n";
	
	s << "S: \n" << ev.S << "\n";
	s << "T: \n" << ev.T << "\n";
	s << "A: \n" << ev.A << "\n";
	s << "B: \n" << ev.B << "\n";	
	s << "nstatev: \n" << ev.nstatev << "\n";
	if (ev.nstatev) {
		s << "statev: \n";
        s << ev.statev.t();
		s << "\n";
	}

	s << "\n\n";

	return s;
}

} //namespace smart
