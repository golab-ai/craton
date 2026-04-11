import numpy
from setuptools import setup, Extension, find_packages
from pathlib import Path
include_dirs = [numpy.get_include()]


def find_files(directory, suffix):
    p = Path(directory)
    return list(map(str, p.glob(f"*.{suffix}")))


libdcd = Extension("compuchem.chemistry.software.gmx_traj_parser.libdcd",
                    ['compuchem/chemistry/software/gmx_traj_parser/libdcd.pyx'],
                    include_dirs = include_dirs + ['compuchem/chemistry/software/gmx_traj_parser/include',
                    'compuchem/chemistry/software/gmx_traj_parser/'])

libxdr = Extension("compuchem.chemistry.software.gmx_traj_parser.libxdr",
                   ['compuchem/chemistry/software/gmx_traj_parser/libxdr.pyx',
                    'compuchem/chemistry/software/gmx_traj_parser/src/xdrfile.c',
                    'compuchem/chemistry/software/gmx_traj_parser/src/xdrfile_xtc.c',
                    'compuchem/chemistry/software/gmx_traj_parser/src/xdrfile_trr.c',
                    'compuchem/chemistry/software/gmx_traj_parser/src/trr_seek.c',
                    'compuchem/chemistry/software/gmx_traj_parser/src/xtc_seek.c',],
                    include_dirs = include_dirs + ['compuchem/chemistry/software/gmx_traj_parser/include',
                    'compuchem/chemistry/software/gmx_traj_parser/'])

lib_align_files = ['compuchem/md_simulation/align/ls_align.pyx']
libalign = Extension("compuchem.md_simulation.align.libalign", lib_align_files, include_dirs=['compuchem/md_simulation/align/LSAlign'],
                     language='c++')

exts = [libdcd, libxdr, libalign]

setup(
    name="compuchem",
    py_modules=["compuchem"],
    # packages = find_packages(),
    version="0.0.1",
    install_requires=[],
    include_package_data = True,
    ext_modules=exts,
    zip_safe = False,
    entry_points={
        'console_scripts': [
            'cc_md = compuchem.command_line.simulation:simulation_app',
            'cc_tool = compuchem.command_line.tools:tools_app',
            'cc_debug = compuchem.command_line.debug:debug_app',
            'cc_benchmark = compuchem.command_line.benchmark:benchmark_app'
        ]
    }
)
