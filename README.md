# SO3: Fast Wigner transforms
[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: http://astro-informatics.github.io/so3/

[![][docs-img]][docs-url]
![Bintray](https://img.shields.io/bintray/v/mdavezac/AstroFizz/so3:AstroFizz?label=bintray%20-%20C%20package)

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
