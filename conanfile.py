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
    requires = "ssht/1.3.2@AstroFizz/stable"
    generators = "cmake"
    exports_sources = [
        "src/c/*",
        "include/*",
        "CMakeLists.txt",
        "cmake/*.cmake",
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
            self._cmake.configure(source_folder=".")
        return self._cmake

    def build(self):
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["so3"]
