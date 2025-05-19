import os
import shutil
from setuptools import setup
from cmake_setuptools import CMakeExtension, CMakeBuildExt

class CustomCMakeBuild(CMakeBuildExt):
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
