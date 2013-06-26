// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details

/*! 
 * \file so3_about.c
 * Print information about the SO3 package, including version
 * and build numbers.
 *
 * Usage: so3_about
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>

int main(int argc, char *argv[]) {

  printf("%s\n", "===================================================================");
  printf("%s\n", "SO3 package to perform Wigner transform on the rotation group SO(3)");
  printf("%s\n", "By Jason McEwen and Yves Wiaux");

  printf("%s\n", "See www.jasonmcewen.org for more information.");
  printf("%s\n", "See LICENSE.txt for license details.");

  printf("%s%s\n", "Version: ", SO3_VERSION);
  printf("%s%s\n", "Build: ", SO3_BUILD);
  printf("%s\n", "===================================================================");

  return 0;

}

