# coding: utf-8

from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy

install_requires = [
    "numpy",
    "cython",
    "biopython",
    "click",
    "mmh3"
]

setup(
    name = "PlasmidPicker",
    version = "0.0.1",
    author = "Masafumi Harada",
    url = "https://github.com/haradama/PlasmidPicker",
    description = "This tool was developed to detect plasmid sequence data in the metagenome.",
    cmdclass = {"build_ext": build_ext},
    ext_modules = cythonize("PlasmidPicker/core.pyx"),
    entry_points = {
        "console_scripts":[
            "plasmidpicker = plasmidpicker.__main__:main",
        ],
    },
    packages =["plasmidpicker"],
    install_requires = install_requires,
    include_dirs = [numpy.get_include()]
)
