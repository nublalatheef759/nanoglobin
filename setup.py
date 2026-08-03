from __future__ import annotations

import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup


extensions = [
    Extension(
        "nanoglobin._genotype_kernel",
        ["nanoglobin/_genotype_kernel.pyx"],
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    )
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "initializedcheck": False,
            "cdivision": True,
        },
    )
)
