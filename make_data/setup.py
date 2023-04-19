import pybind11
from distutils.core import setup, Extension

ext_modules = [
    Extension(
        'CannonBallisticFunctions',
        ['cpp_src/main.cpp'],
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-std=c++20'],
    ),
]

setup(
    name='CannonBallisticFunctions',
    version='0.0.1',
    author='SpaceEye',
    # description='',
    ext_modules=ext_modules,
    requires=['pybind11']
)