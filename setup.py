from skbuild import setup

cmake_args = [
    "-Dpython:BOOL=ON",
    "-Dtests:BOOL=OFF",
]

setup(
    name="pyso3",
    version="1.2.1",
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
            "black",
            "pytest",
        ]
    },
    description="Fast spin spherical transforms",
    url="http://astro-informatics.github.io/ssht/",
    package_dir={"pyso3": "src/pyso3"},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-3",
    packages=["pyso3"],
)
