#ifndef CEDAR_2D_RELAX_GPU_H
#define CEDAR_2D_RELAX_GPU_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>

#include <cedar/2d/gpu/types.h>
#include <cedar/kernels/point_relax.h>
#include <cedar/kernels/line_relax.h>
#include <cedar/2d/gpu/kernel_manager.h>

using real_t = cedar::real_t;
using len_t = cedar::len_t;
#include <src/2d/ftn/mpi/BMG2_SymStd_SETUP_recip.f90.hpp>
#include <src/2d/ftn/mpi/BMG2_SymStd_relax_GS.f90.hpp>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_relax_GS(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                              len_t II, len_t JJ, int kf, int ifd, int nstncl,
	                              int irelax_sym,
	                              int updown, len_t iGs, len_t jGs, void *halof);
	void MPI_BMG2_SymStd_relax_lines_x(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
	void MPI_BMG2_SymStd_relax_lines_y(int k, real_t *SO, real_t *QF, real_t *Q, real_t *SOR,
	                                   real_t *B, len_t II, len_t JJ, len_t iGs, len_t jGs,
	                                   int nog, int nstencil, int irelax_sym, int updown,
	                                   len_t *datadist,
	                                   real_t *rwork, len_t nmsgr, int mpicomm,
	                                   int xlinecomm, int ylinecomm, void *halof);
}

namespace cedar { namespace cdr2 { namespace gpu { namespace mpi {

class rbgs : public kernels::point_relax<stypes>
{
	void setup(const stencil_op<five_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void setup(const stencil_op<nine_pt> & so,
	           relax_stencil & sor) override
	{
		this->setup_impl(so, sor);
	}
	void run(const stencil_op<five_pt> & so,
	         grid_func & x,
	         const grid_func & b,
	         const relax_stencil & sor,
	         cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}
	void run(const stencil_op<nine_pt> & so,
	         grid_func & x,
	         const grid_func & b,
	         const relax_stencil & sor,
	         cycle::Dir cdir) override
	{
		this->run_impl(so, x, b, sor, cdir);
	}

	template<class sten>
	void setup_impl(const stencil_op<sten> & so,
	                relax_stencil & sor)
	{
		int nstencil = stencil_ndirs<sten>::value;
		auto & sod = const_cast<stencil_op<sten>&>(so);
		MPI_BMG2_SymStd_SETUP_recip<ftl::device::GPU>(sod, sor, so.len(0), so.len(1), nstencil);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              grid_func & x,
	              const grid_func & b,
	              const relax_stencil & sor,
	              cycle::Dir cdir)
	{
		int k, kf, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		kf = topo.nlevel();

		if (std::is_same<five_pt, sten>::value)
			ifd = 1;
		else
			ifd = 0;
		nstencil = stencil_ndirs<sten>::value;

		if (cdir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		// ibc = BMG_BCs_definite;

                void* halof = services->fortran_handle<halo_exchange>();

                auto sob = sod.to_buffer();
                auto sorb = sord.to_buffer();
                auto bb = b.to_buffer();
                auto xb = x.to_buffer();

                std::cerr << " == Relaxation == " << std::endl;
                std::cerr << "RHS: " << std::endl << bb << std::endl;
                std::cerr << "Pre-relaxation solution: " << std::endl << xb << std::endl;

                sod.ensure_cpu();
                x.ensure_cpu();
                bd.ensure_cpu();
                sord.ensure_cpu();

		MPI_BMG2_SymStd_relax_GS<ftl::device::GPU>(
                    k, sod, bd, x, sord, so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM, updown, topo.is(0), topo.is(1), halof);

		// MPI_BMG2_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		//                          so.len(0), so.len(1), kf, ifd, nstencil, BMG_RELAX_SYM,
		//                          updown, topo.is(0), topo.is(1),
		//                          this->services->fortran_handle<halo_exchange>());

                // x.mark_cpu_dirty(true);

                // std::cerr << "solution" << std::endl;
                // std::cerr << "has cpu: " << x.has_cpu() << std::endl;
                // std::cerr << "has gpu: " << x.has_gpu() << std::endl;
                // std::cerr << "cpu ptr: " << x.to_flat_buffer().get_host_impl()->get_host_pointer() << std::endl;
                // std::cerr << "dev ptr: " << x.to_flat_buffer().get_dev_impl().get() << std::endl;

                std::cerr << "Post-relaxation solution: " << std::endl << xb << std::endl;
                std::cerr << " ================ " << std::endl;
	}
};

}}}}
#endif
