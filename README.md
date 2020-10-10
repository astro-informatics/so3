# SO3: Fast Wigner transforms
[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: http://astro-informatics.github.io/so3/
[bintray-img]: https://img.shields.io/bintray/v/astro-informatics/astro-informatics/so3:astro-informatics?label=C%20package
[bintray-url]: https://bintray.com/astro-informatics/astro-informatics/so3:astro-informatics/1.3.0:stable/link
[codefactor-img]: https://www.codefactor.io/repository/github/astro-informatics/so3/badge/main
[codefactor-url]: https://www.codefactor.io/repository/github/astro-informatics/so3/overview/main

[![][docs-img]][docs-url]
[![][bintray-img]][bintray-url]
[![][codefactor-img]][codefactor-url]
![CMake Build](https://github.com/astro-informatics/so3/workflows/CMake%20Build/badge.svg)

## DESCRIPTION

The SO3 code provides functionality to perform fast and exact Wigner transforms
using the spherical harmonic transforms from
[SSHT](https://www.github.com/astro-informatics/ssht).

## REFERENCES

 - J. D. McEwen, M. Büttner, B. Leistedt, H. V. Peiris, Y. Wiaux
    [A novel sampling theorem on the rotation group](https://doi.org/10.1109/LSP.2015.2490676)
    IEEE Signal Processing Letters. Vol. 22, No. 12, 2015, pp 2425-2429

## INSTALLATION
The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html):

Both can be installed using pip:

```bash
pip install conan cmake
```

Then so3 can be compiled with:

```bash
git clone http://astro-informatics.github.io/so3/ -b main
mkdir so3/build && cd so3/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -Dconan_deps=True ..
make
make install
```

The above will also download [FFTW](https://www.fftw.org) and
[SSHT](https://www.github.com/astro-informatics/ssht) and compile them, if
necessary.

Instructions for installing the fortran package can be found in docs/index.html.
