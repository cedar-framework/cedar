from spack.package import *

class Cedar(CMakePackage):
    '''The Cedar Framework is a robust, variational multigrid library
    implementing BoxMG on large scale parallel systems.
    '''

    homepage = 'https://github.com/cedar-framework/cedar'
    git = 'https://github.com/cedar-framework/cedar.git'

    maintainers('andrewreisner')

    version('develop', branch='main')

    variant('mpi', default=True, description='Build distributed-memory parallel solver')

    depends_on('blas')
    depends_on('lapack')
    depends_on('boost cxxstd=11 +program_options')
    depends_on('mpi', when='+mpi')

    def cmake_args(self):
        spec = self.spec

        options = [
            self.define_from_variant('ENABLE_MPI', 'mpi')
        ]

        return options
