/*
 * global.hpp
 *
 *  Created on: Feb 17, 2014
 *      Author: Marco
 */

#ifndef GLOBAL_HPP_
#define GLOBAL_HPP_

#include "preprocessor_defines.dat"
#include <math.h>
#include <stdint.h>

template <int dim>
class global;

#include "global_2D.tpp"
#include "global_3D.tpp"

extern const global<2> global2D;
extern const global<3> global3D;


#endif /* GLOBAL_HPP_ */
