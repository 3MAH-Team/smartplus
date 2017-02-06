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

///@file block.cpp
///@brief object that defines a block
///@version 1.0

#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <armadillo>
#include <memory>
#include <smartplus/Libraries/Solver/block.hpp>
#include <smartplus/Libraries/Solver/step_meca.hpp>
#include <smartplus/Libraries/Solver/step_thermomeca.hpp>

using namespace std;
using namespace arma;

namespace smart{

//=====Private methods for ellipsoid_characteristics===================================

//=====Public methods for ellipsoid_characteristics============================================


//@brief default constructor
//-------------------------------------------------------------
block::block()
//-------------------------------------------------------------
{
	number=0;
	nstep=0;
	ncycle=0;
    type=0;
}

/*!
  \brief Constructor
  \param n : number of the block
  \param m : number of steps
  \param k : number of cycles
  \param t : : loading type of steps inside the block
*/

//-------------------------------------------------------------
block::block(int n, int m, int k, int t)
//-------------------------------------------------------------
{
	assert(n>=0);
	assert(m>0);
	assert(k>0);

	number = n;
	nstep = m;
	ncycle = k;
    type = t;
}

/*!
  \brief Constructor with parameters
  \param mnumber : number of the block
  \param mnstep : number of steps
  \param mncycle : number of cycles
  \param mtype : loading type of steps inside the block
  msteps : vector of steps (already defined)
 */

//-------------------------------------------------------------
block::block(int mnumber, int mnstep, int mncycle, int mtype, const vector<shared_ptr<step> > &msteps)
//-------------------------------------------------------------
{	
	assert(mnstep>0);
	assert(mncycle>0);

	number = mnumber;
	nstep = mnstep;
	ncycle = mncycle;
    type = mtype;
    
    steps = msteps;
}

/*!
  \brief Copy constructor
  \param bl block object to duplicate
*/

//------------------------------------------------------
block::block(const block &bl)
//------------------------------------------------------
{  
	number = bl.number;
	nstep = bl.nstep;
	ncycle = bl.ncycle;
    type = bl.type;
    
//    generate();
    
    steps = bl.steps;
}

/*!
  \brief destructor
*/

block::~block() {}

//-------------------------------------------------------------
void block::generate()
//-------------------------------------------------------------
{
	assert(number>=0);
	assert(nstep>0);
	assert(ncycle>0);
	assert(type>0);
    
    switch (type) {
        case 1: {
            
            for (int i=0; i<nstep; i++) {
                shared_ptr<step_meca> sptr_meca(new step_meca);
                steps.push_back(sptr_meca);
            }
            break;
        }
        case 2: {
            
            for (int i=0; i<nstep; i++) {
                shared_ptr<step_thermomeca> sptr_thermomeca(new step_thermomeca);
                steps.push_back(sptr_thermomeca);
            }
            break;
        }
        default: {
            cout << "Please provide a consistent loading type for the block " << number << "\n";
            break;
        }
    }
    
}

/*//-------------------------------------------------------------
void block::initialize()
//-------------------------------------------------------------
{
	assert(number>=0);
	assert(nstep>0);
	assert(ncycle>0);
	assert(type>0);
    
    switch (type) {
        case 1: {
            
            shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(steps[0]);
            sptr_meca->generate(sptr_meca);
            break;
        }
        //        case 2: {
        //
        //break;
        //}
        default: {
            cout << "Please provide a consistent loading type for the block " << number << "\n";
            break;
        }
    }
    
}*/

/*!
  \brief Standard operator = for block
*/

//----------------------------------------------------------------------
block& block::operator = (const block &bl)
//----------------------------------------------------------------------
{
	assert(bl.nstep>0);
	assert(bl.ncycle>0);
  
	number = bl.number;
	nstep = bl.nstep;
	ncycle = bl.ncycle;
    type = bl.type;
    
//    generate();
    
    steps = bl.steps;
    
	return *this;
}

//--------------------------------------------------------------------------
ostream& operator << (ostream &s, const block &bl)
//--------------------------------------------------------------------------
{
	s << "Display info on the block " << bl.number << "\n";
	s << "Number of steps: " << bl.nstep << "\n";
	s << "Number of cycles: " << bl.ncycle << "\n";
	s << "Type of block: " << bl.type << "\n";
	
    switch (bl.type) {
        case 1: {
            
            for (int j=0; j<bl.nstep; j++) {
                shared_ptr<step_meca> sptr_meca = std::dynamic_pointer_cast<step_meca>(bl.steps[j]);
                cout << "\n" << *sptr_meca;
            }
            break;
        }
        case 2: {
            
            for (int j=0; j<bl.nstep; j++) {
                shared_ptr<step_thermomeca> sptr_thermomeca = std::dynamic_pointer_cast<step_thermomeca>(bl.steps[j]);
                cout << "\n" << *sptr_thermomeca;
            }
            break;
        }
        default: {
            cout << "Error Could not find the type of the block \n";
            break;
        }
    }
	return s;
}

} //namespace smart
