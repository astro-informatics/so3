from conans import CMake, ConanFile


class So3Conan(ConanFile):
    name = "so3"
    version = "1.3.0"
    license = "GPL-3"
    url = "https://github.com/astro-informatics/so3"
    homepage = "https://github.com/astro-informatics/so3"
    description = "Fast and accurate Wigner transforms"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"fPIC": [True, False]}
    default_options = {"fPIC": True}
    requires = "ssht/1.3.2@astro-informatics/stable"
    generators = "cmake"
    exports_sources = [
        "src/c/*",
        "include/*",
        "CMakeLists.txt",
        "cmake/*.cmake",
        "tests/*.c",
        "tests/*.h",
        "tests/CMakeLists.txt",
    ]

    def configure(self):
        if self.settings.compiler == "Visual Studio":
            del self.options.fPIC
        self.options["ssht"].fPIC = self.options.fPIC
        del self.settings.compiler.libcxx

    @property
    def cmake(self):
        if not hasattr(self, "_cmake"):
            self._cmake = CMake(self)
            self._cmake.definitions["tests"] = True
            self._cmake.definitions["conan_deps"] = True
            self._cmake.definitions["fPIC"] = self.options.fPIC
            self._cmake.configure(build_folder="build")
        return self._cmake

    def build(self):
        from pathlib import Path

        path = Path(self.source_folder)
        build = Path(self.source_folder) / "build"
        build.mkdir(exist_ok=True)
        (path / "conanbuildinfo.cmake").rename(path / "build" / "conanbuildinfo.cmake")
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["so3"]
