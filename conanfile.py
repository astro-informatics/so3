from conans import ConanFile, CMake


class So3Conan(ConanFile):
    name = "so3"
    version = "1.2.1"
    license = "GPL-3"
    url = "https://github.com/astro-informatics/so3"
    homepage = "https://github.com/astro-informatics/so3"
    description = "Fast and accurate Wigner transforms"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"fPIC": [True, False]}
    default_options = {"fPIC": True}
    requires = "ssht/1.3.1@AstroFizz/stable"
    generators = "cmake"
    exports_sources = [
        "src/c/*",
        "include/*",
        "CMakeLists.txt",
        "cmake/*.cmake",
    ]

    def configured_cmake(self):
        cmake = CMake(self)
        cmake.definitions["tests"] = True
        cmake.definitions["conan_deps"] = True
        cmake.definitions["fPIC"] = self.options.fPIC
        cmake.configure(source_folder=".")
        return cmake

    def build(self):
        cmake = self.configured_cmake()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = self.configured_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["so3"]
