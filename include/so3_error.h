// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*! \file so3_error.h
 *  Error macros used in SO3 package.
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#ifndef SO3_ERROR
#define SO3_ERROR

#include <stdio.h>

// Put this macro in a block so that it can be used with single-line
// if-statements.
#define SO3_ERROR_GENERIC(comment) 					\
{                                                                       \
  printf("ERROR: %s.\n", comment);					\
  printf("ERROR: %s <%s> %s %s %s %d.\n",				\
	 "Occurred in function",					\
	   __PRETTY_FUNCTION__,						\
	   "of file", __FILE__,						\
	   "on line", __LINE__);					\
  exit(1);                                                              \
}

#define SO3_ERROR_MEM_ALLOC_CHECK(pointer)				\
  if(pointer == NULL) {							\
    SO3_ERROR_GENERIC("Memory allocation failed")			\
  }

#endif
