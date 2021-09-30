from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os

# If the user set GSL_PREFIX, use it. Otherwise check some standard locations.
gsl_prefix = ""
try:
    gsl_prefix = os.environ['GSL_PREFIX']
except:
    print("Environment variable GSL_PREFIX not set. Trying known locations...")
    prefixes = ["/usr/", "/usr/local/", "/opt/local/", "/sw/"]

    for path in prefixes:
        print("Checking '%s'..." % path)
        try:
            os.stat(path + "include/gsl/gsl_odeiv2.h")
            gsl_prefix = path
            print("Found GSL in '%s'" % gsl_prefix)
            break
        except:
            pass

if gsl_prefix == "":
    print("Could not find GSL. Stopping...")
    import sys

    sys.exit(1)


# If the user set NO_OPENMP, proceed with these options. Otherwise add options clang uses.
libraries=['gsl', 'gslcblas']
library_dirs=[gsl_prefix + "/lib"]
extra_compile_args=["-O3", "-ffast-math", "-Wall"]
extra_link_args=[]
try:
    os.environ["NO_OPENMP"]
except:
    extra_compile_args.append('-fopenmp')
    libraries.append('gomp')
    # library_dirs.append("/opt/local/lib/libomp")

# Define the extension
extension = Extension("pism_dbg",
                      sources=["python/pism_dbg.pyx",
                               "src/upslope_area.cc",
                               "src/accumulated_flow.cc",
                               "src/initialize_mask.cc",
                               "src/DEM.cc"
                               ],
                      include_dirs=[numpy.get_include(), 'src', gsl_prefix + '/include'],
                      library_dirs=library_dirs,
                      libraries=libraries,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args,
                      language="c++")

setup(
    name = "pism_dbg",                       # "drainage basin generator"
    version = "0.2.0",
    description = "PISM drainage basin generator",
    long_description = """
    pism_dbg is the 'drainage basin generator' for regional modeling using PISM.
    See http://www.pism-docs.org for details.""",
    author = "PISM authors",
    author_email = "uaf-pism@alaska.edu",
    url = "https://github.com/pism/regional-tools",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Environment :: X11 Applications",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Topic :: Utilities",
    ],
    cmdclass={"build_ext": build_ext},
    ext_modules=[extension],
)
