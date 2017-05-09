# Cedar

The Cedar Framework is a robust, variational multigrid library implementing
BoxMG on large scale parallel systems.

The Cedar Framework is a collaborative project between the University of
Illinois at Urbana-Champaign and Los Alamos National Laboratory.  Work on this
project is based in part upon work supported by the Department of Energy,
National Nuclear Security Administration, under Award Number DE-NA0002374.

TODO logo

Cedar is release under the BSD-3-Clause (see LICENSE.txt).

# Installation

```
mkdir build
cd build
ccmake ..
```

1. `[c]` to configure.  As a basic build, use
    - `CMAKE_INSTALL_PREFIX` as `./local_install/`
    - `MPI_LIBRARY` is set correctly

2. `[c]` to configure again
3. `[g]` to generate
4. `make` to compile (or `make -j4` to compile with 4 threads)
5. `make install` 

## CMAKE Options

ENABLE_3D           | Build with 3D support
ENABLE_EXAMPLES     | Build examples
ENABLE_MPI          | Build with MPI support
ENABLE_UNIT_TESTS   | Build unit tests

# Examples

Following the installation instructions above, run the following from the build directory:
- `./examples/`

# Citing

<pre>
@misc{cedar,
  author = {David Moulton and Luke N. Olson and Andrew Reisner},
  title = {Cedar Framework},
  year = {2017},
  url = {https://github.com/cedar-framework/cedar},
  note = {Version 0.1},
}
</pre>

# Getting Help

Creat an [issue](https://github.com/cedar-framework/cedar/issues).
