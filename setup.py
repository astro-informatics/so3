from skbuild import setup

cmake_args = [
    "-Dpython:BOOL=ON",
    "-Dtests:BOOL=OFF",
]

setup(
    name="so3",
    version="1.3.0",
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
    package_dir={"so3": "src/so3"},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-3",
    packages=["so3"],
)
