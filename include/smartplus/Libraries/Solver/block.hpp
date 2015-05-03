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

///@file block.hpp
///@brief object that defines a block
///@version 1.0

#pragma once
#include <iostream>
#include <armadillo>
#include "step.hpp"
#include <memory>

using namespace std;
using namespace arma;

namespace smart{

//======================================
class block
//======================================
{
	private:

	protected:

	public :
		int number;
        int nstep;
        int ncycle;
        int type;
    
        vector<shared_ptr<step> > steps;
    
        block(); 	//default constructor
		block(int, int, int, int);	//constructor - allocates memory for the step vector
        block(int, int, int, int, const std::vector<shared_ptr<step> > &); //Constructor with parameters
		block(const block&);	//Copy constructor
		~block();
		
		void generate();
//		void initialize();
    
		virtual block& operator = (const block&);
		
        friend ostream& operator << (ostream&, const block&);
};

} //namespace smart
