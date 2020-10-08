// S03 package to perform Wigner transform on the rotation group SO(3)
// Copyright (C) 2013 Martin Büttner and Jason McEwen
// See LICENSE.txt for license details

/*!
 * \file so3_about.c
 * Print information about the SO3 package, including version
 * and build numbers.
 *
 * Usage: so3_about
 *
 * \author <a href="mailto:m.buettner.d@gmail.com">Martin Büttner</a>
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>
#ifdef BUILT_WITH_CMAKE
#include "so3_version.h"
#endif

int main(int argc, char *argv[]) {

#ifdef BUILT_WITH_CMAKE
  printf("%s", so3_info());
#else
  printf("%s\n", "===================================================================");
  printf("%s\n", "SO3 package to perform Wigner transform on the rotation group SO(3)");
  printf("%s\n", "By Martin Büttner and Jason McEwen");

  printf("%s\n", "See www.jasonmcewen.org for more information.");
  printf("%s\n", "See LICENSE.txt for license details.");

  printf("%s%s\n", "Version: ", SO3_VERSION);
  printf("%s%s\n", "Build: ", SO3_BUILD);
  printf("%s\n", "===================================================================");
#endif

  return 0;
}
