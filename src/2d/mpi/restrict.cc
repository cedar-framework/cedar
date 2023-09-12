#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/2d/mpi/restrict.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/mpi/BMG2_SymStd_restrict.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_restrict(int kf, int kc, int nog,
	                              real_t *Q, real_t *QC, real_t *CI,
	                              len_t nx, len_t ny, len_t nxc, len_t nyc,
	                              len_t iGs, len_t jGs);
}


template <typename T>
static void print_cuda_buffer(ftl::Buffer<T> buf) {
    const long numel = buf.get_numel();
    const long bytes = sizeof(T) * numel;
    auto cuda_buffer =
        std::dynamic_pointer_cast<ftl::CUDABaseBuffer<double, unsigned int>>(
            buf.get_dev_impl());
    assert(cuda_buffer != nullptr);
    double* device_pointer = cuda_buffer->get_device_pointer();
    std::unique_ptr<double[]> host_buffer = std::make_unique<double[]>(numel);

    cudaMemcpy(host_buffer.get(), device_pointer, bytes, cudaMemcpyDeviceToHost);

    for (std::size_t i = 0; i < numel; ++i) {
        std::cerr << host_buffer[i] << " ";
    }
    std::cerr << std::endl;
}

template <typename T>
static void print_host_buffer(ftl::Buffer<T> buf) {
    const long numel = buf.get_numel();
    const long bytes = sizeof(T) * numel;
    auto host_buffer =
        std::dynamic_pointer_cast<ftl::BaseHostBufferType<double, unsigned int>>(
            buf.get_host_impl());
    assert(host_buffer != nullptr);
    double* host_pointer = host_buffer->get_host_pointer();

    for (std::size_t i = 0; i < numel; ++i) {
        std::cerr << host_pointer[i] << " ";
    }
    std::cerr << std::endl;
}

namespace cedar { namespace cdr2 { namespace mpi {

	void restrict_f90::update_periodic(mpi::grid_func & q,
	                            const grid_topo & topof,
	                            const grid_topo & topoc,
	                            std::array<bool, 3> periodic)
	{
		std::array<len_t, 2> ngf({topof.nglobal(0), topof.nglobal(1)});
		std::array<len_t, 2> ngc({topoc.nglobal(0), topoc.nglobal(1)});

                q.ensure_cpu();

		if (periodic[0] and ((ngf[0]/2 + 1) != ngc[0])) {
			if (topof.coord(0) == 0) {
				for (auto j : q.grange(1)) {
					q(0,j) = 0;
				}
			}
			if (topof.coord(0) == (topof.nproc(0) - 1)) {
				for (auto j : q.grange(1)) {
					q(q.len(0)-1,j) = 0;
				}
			}
		}

		if (periodic[1] and ((ngf[1]/2 + 1) != ngc[1])) {
			if (topof.coord(1) == 0) {
				for (auto i : q.grange(0)) {
					q(i,0) = 0;
				}
			}
			if (topof.coord(1) == (topof.nproc(1) - 1)) {
				for (auto i : q.grange(0)) {
					q(i,q.len(1)-1) = 0;
				}
			}
		}

                q.mark_cpu_dirty(true);
	}


	void restrict_f90::run(const restrict_op & R,
                       const grid_func & fine,
                       grid_func & coarse)
	{
		int kf, kc, nog;
		auto & Rd = const_cast<restrict_op&>(R);
		prolong_op & P = Rd.getP();
		grid_topo & topo = P.grid();
		const grid_topo & fine_topo = fine.grid();
		const grid_topo & coarse_topo = coarse.grid();
		auto & fined = const_cast<mpi::grid_func&>(fine);

		nog = kf = topo.nlevel();
		kc = topo.level();

		// conditionally zero periodic entries
		update_periodic(fined, fine_topo, coarse_topo, params->periodic);

                if (fined.has_gpu() || coarse.has_gpu() || P.has_gpu()) {
                    fined.mark_cpu_dirty(true);
                    fined.ensure_gpu();
                    coarse.ensure_gpu();
                    P.ensure_gpu();

                    MPI_BMG2_SymStd_restrict<ftl::device::GPU>(
                        kf, kc, nog, fined, coarse, P,
                        fined.len(0), fined.len(1),
                        coarse.len(0), coarse.len(1),
                        fine_topo.is(0), fine_topo.is(1));

                    coarse.ensure_cpu();
                    std::cerr << "GPU: " << std::endl << coarse << std::endl;
                } //else {
                {
                    fined.ensure_cpu();
                    coarse.ensure_cpu();
                    P.ensure_cpu();

                    MPI_BMG2_SymStd_restrict(kf, kc, nog,
                                             fined.data(), coarse.data(), P.data(),
                                             fined.len(0), fined.len(1),
                                             coarse.len(0), coarse.len(1),
                                             fine_topo.is(0), fine_topo.is(1));

                    coarse.mark_cpu_dirty(true);

                    std::cerr << "CPU: " << std::endl << coarse << std::endl;
                }
	}

}}}
