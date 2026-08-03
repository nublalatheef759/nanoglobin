from __future__ import annotations

from pathlib import Path

from setuptools import Extension, setup


ROOT = Path(__file__).resolve().parent


def build_extensions():
    """Build the typed genotype kernel while keeping a Python fallback."""

    try:
        import numpy
        from Cython.Build import cythonize
    except ImportError:
        # PEP 517 installs these through pyproject.toml.  Returning no extension
        # keeps source-tree tooling usable in deliberately minimal environments;
        # production/editable installs build the Cython path.
        return []

    extensions = [
        Extension(
            "nanoglobin._genotype_fast",
            [str(ROOT / "nanoglobin" / "_genotype_fast.pyx")],
            include_dirs=[numpy.get_include()],
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
