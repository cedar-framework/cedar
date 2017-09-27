#ifndef CEDAR_3D_RELAX_MPI_H
#define CEDAR_3D_RELAX_MPI_H

#include <type_traits>

#include "cedar/2d/ftn/BMG_parameters_c.h"
#include <cedar/kernel_params.h>
#include <cedar/halo_exchanger.h>
#include <cedar/cycle/types.h>
#include <cedar/3d/mpi/grid_func.h>
#include <cedar/3d/mpi/stencil_op.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_relax_GS(int kg, real_t *so, real_t *qf, real_t *q, real_t *sor,
	                              len_t nlx, len_t nly, len_t nlz, len_t ngx, len_t ngy, len_t ngz,
	                              int nog, int ifd, int nstencil, int nsorv,
	                              int irelax_sym, int updown,
	                              len_t igs, len_t jgs, len_t kgs, void *halof);
}


namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	template<class sten>
	void mpi_relax_rbgs_point(const kernel_params & params,
	                          halo_exchanger_base *halof,
	                          const mpi::stencil_op<sten> & so,
	                          mpi::grid_func & x,
	                          const mpi::grid_func & b,
	                          const relax_stencil & sor,
	                          cycle::Dir cycle_dir)
	{
		using namespace cedar::cdr3;
		int k, ifd;
		int updown, nstencil;

		auto & sod = const_cast<mpi::stencil_op<sten>&>(so);
		grid_topo & topo = sod.grid();
		relax_stencil & sord = const_cast<relax_stencil&>(sor);
		mpi::grid_func & bd = const_cast<mpi::grid_func&>(b);

		k = topo.level()+1;
		nstencil = stencil_ndirs<sten>::value;

		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		if (cycle_dir == cycle::Dir::UP) updown = BMG_UP;
		else updown = BMG_DOWN;

		int nsorv = 2;

		// ibc = BMG_BCs_definite;
		MPI_BMG3_SymStd_relax_GS(k, sod.data(), bd.data(), x.data(), sord.data(),
		                         so.len(0), so.len(1), so.len(2),
		                         topo.nglobal(0), topo.nglobal(1), topo.nglobal(2),
		                         topo.nlevel(), ifd, nstencil, nsorv,
		                         BMG_RELAX_SYM, updown,
		                         topo.is(0), topo.is(1), topo.is(2),
		                         halof);
	}
}

}}}

#endif
