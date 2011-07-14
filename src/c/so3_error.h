// SO3 package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


/*! \file so3_error.h
 *  Error macros used in SO3 package.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#ifndef SO3_ERROR
#define SO3_ERROR


#include <stdio.h>


#define SO3_ERROR_GENERIC(comment) 					\
  printf("ERROR: %s.\n", comment);					\
  printf("ERROR: %s <%s> %s %s %s %d.\n",				\
	 "Occurred in function",					\
	   __PRETTY_FUNCTION__,						\
	   "of file", __FILE__,						\
	   "on line", __LINE__);					\
  exit(1);

#define SO3_ERROR_MEM_ALLOC_CHECK(pointer)				\
  if(pointer == NULL) {							\
    SO3_ERROR_GENERIC("Memory allocation failed")			\
  }


#endif
