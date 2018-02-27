#include "cedar/2d/ftn/BMG_parameters_c.h"
#include "cedar/2d/relax/setup_relax.h"

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_recip(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_x(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
	void MPI_BMG2_SymStd_SETUP_lines_y(real_t *SO, real_t *SOR, len_t Nx, len_t Ny, int NStncl);
}


namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	using namespace cedar::cdr2;

	template<>
	void mpi_setup_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op<five_pt> & so,
	                          relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);

		nstencil = 3;

		MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}


	template<>
	void mpi_setup_rbgs_point(const kernel_params & params,
	                          const mpi::stencil_op<nine_pt> & so,
	                          relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);

		nstencil = 5;

		MPI_BMG2_SymStd_SETUP_recip(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}


	template<>
	void mpi_setup_rbgs_x(const kernel_params & params,
	                      const mpi::stencil_op<five_pt> & so,
	                      relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);


		nstencil = 3;

		MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}


	template<>
	void mpi_setup_rbgs_x(const kernel_params & params,
	                      const mpi::stencil_op<nine_pt> & so,
	                      relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);

		nstencil = 5;

		MPI_BMG2_SymStd_SETUP_lines_x(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}


	template<>
	void mpi_setup_rbgs_y(const kernel_params & params,
	                      const mpi::stencil_op<five_pt> & so,
	                      relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<five_pt>&>(so);

		nstencil = 3;

		MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}


	template<>
	void mpi_setup_rbgs_y(const kernel_params & params,
	                      const mpi::stencil_op<nine_pt> & so,
	                      relax_stencil & sor)
	{
		int nstencil;

		auto & sod = const_cast<mpi::stencil_op<nine_pt>&>(so);

		nstencil = 5;

		MPI_BMG2_SymStd_SETUP_lines_y(sod.data(), sor.data(), so.len(0), so.len(1), nstencil);
	}
}

}}}
