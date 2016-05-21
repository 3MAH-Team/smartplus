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

///@file ODF.hpp
///@brief Characteristics of an ODF
///@version 0.9

#pragma once

#include <iostream>
#include <string>
#include <armadillo>
#include "peak.hpp"

namespace smart{

//======================================
class ODF
//======================================
{
	private:

	protected:

	public :

        int Nphases; //Number of phases for the ODF discretization
        int Angle;  //Angle for the ODF
        bool radian;
        int n_densities; //number of angle to evaluate the densities (mostly for plotting purposes)
    
        std::vector<peak> peaks; //peaks
        double norm;
    
        arma::vec densities; //sum of densities
        //mat full_ODF stores the ODF densities for each angle (rows), for each peak (column)
        arma::vec limits; //minimal and maximal angles of orientation
    
		ODF(); 	//default constructor
    
        ODF(const int &, const int &, const int &, const std::vector<peak> &, const double &, const arma::vec &, const arma::vec &, const bool & = false);

		ODF(const ODF&);	//Copy constructor
        virtual ~ODF();

		virtual ODF& operator = (const ODF&);
    
        friend std::ostream& operator << (std::ostream&, const ODF&);
};

} //namespace smart