from skbuild import setup

cmake_args = [
    "-DBUILD_TESTING:BOOL=OFF",
    "-Dconan_deps=ON",
    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
]

setup(
    name="so3",
    version="1.3.4",
    author="Jason McEwen",
    install_requires=["numpy", "scipy"],
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
    description="Fast and exact Wigner Transforms",
    url="http://astro-informatics.github.io/so3/",
    package_dir={"so3": "src/so3"},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-3",
    packages=["so3"],
)
