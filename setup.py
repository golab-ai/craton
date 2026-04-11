# import numpy
from setuptools import setup, Extension, find_packages
from pathlib import Path
# include_dirs = [numpy.get_include()]


def find_files(directory, suffix):
    p = Path(directory)
    return list(map(str, p.glob(f"*.{suffix}")))



setup(
    name="craton",
    py_modules=["craton"],
    descript="the testing package of craton ",
    author="CFL",
    # packages = find_packages(),
    entry_points={
        'console_scripts': [
            'craton = craton.command_line:craton_app',
        ]
    },
    version="0.0.1",
    install_requires=[],
    zip_safe = False,
)
