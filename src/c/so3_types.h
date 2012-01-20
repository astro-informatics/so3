// SO3 package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details

/*! \mainpage SO3 C documentation
 *
 * The SO3 code provides functionality to perform fast and exact spin
 * spherical harmonic transforms based on the sampling theorem on the
 * sphere derived in our paper: <i>A novel sampling theorem on the
 * sphere</i> (<a href="http://arxiv.org/abs/XXX.XXX">ArXiv</a>|
 * <a href="http://dx.doi.org/10.1111/XXX">DOI</a>).
 *   
 * We document the C source code here.  For an example of usage, 
 * see the so3_test.c program.
 * For installation instructions, see the general SO3 
 * documentation available 
 * <a href="../../index_so3.html">here</a>.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 * \version 0.1
 */


/*! \file so3_types.h
 *  Types used in SO3 package.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#ifndef SO3_TYPES
#define SO3_TYPES

#define SO3_PI    3.141592653589793238462643383279502884197
#define SO3_PION2 1.570796326794896619231321691639751442099

#define SO3_SQRT2 1.41421356237309504880168872420969807856967


#define SO3_VERSION 0.1
#define SO3_PROMPT "[so3-0.1]"



#endif
