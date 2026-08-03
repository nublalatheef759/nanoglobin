from __future__ import annotations

from setuptools import Extension, setup


def build_extensions():
    """Build the typed genotype kernel while keeping a Python fallback."""

    try:
        import numpy
        from Cython.Build import cythonize
    except ImportError:
        # PEP 517 installs these through pyproject.toml. Returning no extension
        # keeps source-tree tooling usable in deliberately minimal environments;
        # production/editable installs build the Cython path.
        return []

    extensions = [
        Extension(
            "nanoglobin._genotype_fast",
            ["nanoglobin/_genotype_fast.pyx"],
            include_dirs=[numpy.get_include()],
            define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        )
    ]
    return cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "boundscheck": False,
            "wraparound": False,
            "initializedcheck": False,
            "cdivision": True,
        },
    )


setup(ext_modules=build_extensions())
