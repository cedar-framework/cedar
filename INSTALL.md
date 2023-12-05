# Install instructions for GPU acceleration

## Prerequisites

- Make sure you have access to the [GPU translation layer on Github](https://github.com/cedar-framework/fortran-to-loopy)
- Make sure you have at least Python 3.10 installed on your computer.  I have found success using Miniconda3 on clusters to get a more recent Python version installed.

## Install steps

1. If you are using a Python venv, activate it now.  Download [this fork of loopy](https://github.com/nicknytko/loopy)
   and install it using `pip install .` when in the directory.
2. Clone cedar and check out the `loopy-gpu` branch.
3. Make a `build` folder in the checked out source tree and build with CMake as usual with `cmake ..`.  It should automatically pull in dependencies like Tausch and the translation layer.
4. Build with `make`.
5. Once built, the `examples/mpi-poisson-2d` example is set up to use GPU acceleration.  Make sure you set `"gpu": true` in the `config.json` file or pass in the `--gpu` command-line argument to enable GPU support.
