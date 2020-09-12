# SO3: Fast Wigner transforms
[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: http://astro-informatics.github.io/so3/
[bintray-img]: https://img.shields.io/bintray/v/mdavezac/AstroFizz/so3:AstroFizz?label=C%20package
[bintray-url]: https://bintray.com/mdavezac/AstroFizz/so3:AstroFizz/1.3.1:stable/link

[![][docs-img]][docs-url]
[![][bintray-img]][bintray-url]

## DESCRIPTION

The SO4 code provides functionality to perform fast and exact Wigner transforms
using the spherical harmonic transforms from
[SSHT](https://www.github.com/astro-informatics/ssht).

## VERSION
Release 1.2.0, Sept 20

## REFERENCES

 - J. D. McEwen, M. Büttner, B. Leistedt, H. V. Peiris, Y. Wiaux
    [A novel sampling theorem on the rotation group](https://doi.org/10.1109/LSP.2015.2490676)
    IEEE Signal Processing Letters. Vol. 22, No. 12, 2015, pp 2425-2429

## INSTALLATION

The can be installed using [CMake](https://cmake.org):

```cmake
git clone http://astro-informatics.github.io/so3/ -b main
mkdir so3/build && cd so3/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make install
```

The above will also download [FFTW](https://www.fftw.org) and
[SSHT](https://www.github.com/astro-informatics/ssht) and compile them, if
necessary.
