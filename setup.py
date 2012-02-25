from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
import os

# If the user set GSL_PREFIX, use it. Otherwise check some standard locations.
prefix = ""
try:
    prefix = os.environ['GSL_PREFIX']
except:
    print "Environment variable GSL_PREFIX not set. Trying known locations..."
    prefixes = ["/usr/", "/usr/local/", "/opt/local/", "/sw/"]

    for path in prefixes:
        print "Checking '%s'..." % path
        try:
            os.stat(path + "include/gsl/gsl_odeiv.h")
            prefix = path
            print "Found GSL in '%s'" % prefix
            break
        except:
            pass

if prefix == "":
    print "Could not find GSL. Stopping..."
    import sys
    sys.exit(1)


# If the user set NO_OPENMP, proceed with these options. Otherwise add options GCC uses.
libraries=['gsl', 'gslcblas']
extra_compile_args=["-O3", "-ffast-math", "-Wall"]
try:
    os.environ['NO_OPENMP']
except:
    libraries.append('gomp')
    extra_compile_args.append('-fopenmp')

# Define the extension
extension = Extension("dbg",
                      sources=["python/dbg.pyx",
                               "src/upslope_area.cc",
                               "src/accumulated_flow.cc",
                               "src/initialize_mask.cc",
                               "src/DEM.cc"
                               ],
                      include_dirs=[numpy.get_include(), 'src', prefix + '/include'],
                      library_dirs=[prefix + "/lib"],
                      libraries=libraries,
                      extra_compile_args=extra_compile_args,
                      language="c++")

setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [extension]
      )
