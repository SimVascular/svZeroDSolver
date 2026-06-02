import os
import shutil
import subprocess
import sys
from setuptools import setup
from cmake_setuptools import CMakeExtension, CMakeBuildExt

class CustomCMakeBuild(CMakeBuildExt):
    def build_extension(self, ext):
        # cmake_setuptools' default build_extension hard-codes `make`, which
        # does not exist on Windows. Reimplement the configure + build steps
        # with a generator-agnostic `cmake --build` so the same code path
        # works on Linux, macOS, and Windows.
        if not isinstance(ext, CMakeExtension):
            return super().build_extension(ext)

        cmake = os.environ.get("CMAKE_EXE", shutil.which("cmake"))
        if not cmake:
            raise RuntimeError(
                "cmake executable not found. Set CMAKE_EXE or update your PATH"
            )

        output_dir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        build_type = "Debug" if self.debug else "Release"

        configure = [
            cmake,
            ext.sourcedir,
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + output_dir,
            "-DCMAKE_BUILD_TYPE=" + build_type,
        ]
        # On Windows, build with Ninja (a declared build dependency) to match
        # the known-good build recipe; an explicit CMAKE_GENERATOR wins. Also
        # disable Fortran: Eigen runs `enable_language(Fortran OPTIONAL)`, and
        # under Ninja that picks up an incompatible gfortran from the runner's
        # PATH (MinGW/Strawberry) and fails to link. We have no Fortran sources.
        if sys.platform == "win32":
            if not os.environ.get("CMAKE_GENERATOR"):
                configure += ["-G", "Ninja"]
            configure += ["-DCMAKE_Fortran_COMPILER:FILEPATH="]
        configure += [
            x for x in os.environ.get("CMAKE_COMMON_VARIABLES", "").split(" ") if x
        ]

        os.makedirs(self.build_temp, exist_ok=True)
        subprocess.check_call(configure, cwd=self.build_temp)
        subprocess.check_call(
            [cmake, "--build", ".", "--target", ext.name,
             "--config", build_type, "--parallel"],
            cwd=self.build_temp,
        )

    def run(self):
        # -------------------------------------------------
        # 1. Build the C++ extension *without* the default
        #    setuptools copy step (set inplace False)
        # -------------------------------------------------
        inplace_orig = self.inplace        # remember
        self.inplace = False               # inhibit copy_extensions_to_source
        super().run()                      # runs CMake
        self.inplace = inplace_orig        # restore flag

        # -------------------------------------------------
        # 2. Locate the compiled library
        # -------------------------------------------------
        build_temp = os.path.abspath(self.build_temp)
        search_root = os.path.join(build_temp, "python")
        dest_dir = os.path.dirname(self.get_ext_fullpath("pysvzerod"))

        for root, _, files in os.walk(search_root):
            for f in files:
                if f.startswith("pysvzerod") and f.endswith((".so", ".pyd", ".dll", ".dylib")):
                    src = os.path.join(root, f)
                    os.makedirs(dest_dir, exist_ok=True)
                    shutil.copy2(src, os.path.join(dest_dir, f))
                    print(f"[INFO] copied {src} -> {dest_dir}")
                    return

        raise RuntimeError("pysvzerod binary not found in build tree")

setup(
    ext_modules=[CMakeExtension("pysvzerod")],
    cmdclass={"build_ext": CustomCMakeBuild},
)
