from skbuild import setup

cmake_args = [
    "-Dpython:BOOL=ON",
    "-Dtests:BOOL=OFF",
    "-Dconan_fftw=ON",
]

setup(
    name="pyso3",
    version="2.0",
    author="Jason McEwen",
    install_requires=["numpy", "cython", "scipy"],
    extras_require={
        "dev": [
            "setuptools",
            "wheel",
            "scikit-build",
            "cmake",
            "ninja",
            "cython",
            "conan",
        ]
    },
    description="Fast spin spherical transforms",
    url="http://astro-informatics.github.io/ssht/",
    package_dir={"pyso3": "src/pyso3"},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-2",
    packages=["pyso3"],
)
