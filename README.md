# SO3: Fast Wigner transforms
[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: http://astro-informatics.github.io/so3/
[conan-img]: https://img.shields.io/badge/ConanCenter-C%20Package-red.svg
[conan-url]: https://conan.io/center/astro-informatics-so3
[pypi-img]: https://badge.fury.io/py/so3.svg
[pypi-url]: https://badge.fury.io/py/so3
[codefactor-img]: https://www.codefactor.io/repository/github/astro-informatics/so3/badge/main
[codefactor-url]: https://www.codefactor.io/repository/github/astro-informatics/so3/overview/main

[![][docs-img]][docs-url]
[![][conan-img]][conan-url]
[![][pypi-img]][pypi-url]
[![][codefactor-img]][codefactor-url]
![CMake Build](https://github.com/astro-informatics/so3/workflows/CMake%20Build/badge.svg)
![Python Build](https://github.com/astro-informatics/so3/workflows/Python%20Build/badge.svg)

## DESCRIPTION

The <strong>SO3</strong> code provides functionality to perform fast and exact Wigner transforms based on the sampling theorem on the rotation group derived in [McEwen et al. (2015)](http://www.jasonmcewen.org/publication/mcewen-so-3/).

## INSTALLATION

 The python package, <strong>so3</strong> (pyso3 was taken), is available on <a href="https://pypi.org/project/so3/">pypi</a> and can be installed with: 
 
 ```bash
 pip install so3
 ```

The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html):

Both can be installed using pip:

```bash
pip install conan cmake
```

Then **SO3** can be compiled with:

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


## DOCUMENTATION

Further documentation is available [here](https://astro-informatics.github.io/so3/).

Usage for the python package is also given in the package docstring.


## REFERENCING

If you use **SO3** for work that results in publication, please reference <a
href="http://github.com/astro-informatics/so3">https://github.com/astro-informatics/so3/</a>
and cite our related academic paper:

- J. D. McEwen, M. B&uuml;ttner, B. Leistedt, H. V. Peiris, Y. Wiaux, [A novel sampling theorem on the rotation group](http://www.jasonmcewen.org/publication/mcewen-so-3/), IEEE Sig. Proc. Let., 22(12):2425-2429, 2015 [(arXiv:1508.03101)](https://arxiv.org/abs/1508.03101")



## LICENSE

SO3 is released under the GPL-3 license.  For further details see 
[LICENSE](https://github.com/astro-informatics/so3/blob/main/LICENSE).

## AUTHORS

**SO3** was initially developed by Martin B&uuml;ttner, <a href="http://www.jasonmcewen.org/">Jason McEwen</a>, and Boris Leistedt but significant contributors have since been made by a number of [others](https://github.com/astro-informatics/so3/graphs/contributors).
  </p>