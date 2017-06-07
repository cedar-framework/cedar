#ifndef CEDAR_3D_INTER_MPI_SETUP_INTERP_H
#define CEDAR_3D_INTER_MPI_SETUP_INTERP_H

#include <type_traits>

#include <cedar/kernel_params.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/inter/mpi/prolong_op.h>
#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/3d/mpi/halo.h>


extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_interp_OI(int kgf, int kgc, real_t *so, real_t *soc,
	                                     real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     int nog, int ifd, int nstencil, int irelax, real_t *yo,
	                                     int nogm, len_t *IGRD, len_t *iwork, len_t NMSGi,
	                                     int *pMSG, real_t *buffer, len_t nmsgr, int myproc,
	                                     int mpicomm, int jpn);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr3 { namespace kernel {

namespace impls
{
	namespace mpi = cedar::cdr3::mpi;

	template<class sten>
		void store_fine_op(mpi::stencil_op<sten> & fop,
		                   inter::mpi::prolong_op & P);

	template<>
		inline void store_fine_op(mpi::stencil_op<seven_pt> & fop,
		                          inter::mpi::prolong_op & P)
	{
		P.fine_op_seven = &fop;
		P.fine_is_seven = true;
	}

	template<>
		inline void store_fine_op(mpi::stencil_op<xxvii_pt> & fop,
		                          inter::mpi::prolong_op & P)
	{
		P.fine_op_xxvii = &fop;
		P.fine_is_seven = false;
	}

	template<class sten>
	void mpi_setup_interp(const kernel_params & params,
	                      const mpi::stencil_op<sten> & fop,
	                      const mpi::stencil_op<xxvii> & cop,
	                      inter::mpi::prolong_op & P)
	{
		int ifd, nstencil;

		auto & fopd = const_cast<mpi::stencil_op<sten>&>(fop);
		auto & copd = const_cast<mpi::stencil_op<xxvii>&>(cop);
		grid_topo & topo = fopd.grid();
		MsgCtx *ctx = (MsgCtx*) fopd.halo_ctx;

		store_fine_op(fopd, P);

		nstencil = stencil_ndirs<sten>::value;

		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		int rank;
		MPI_Comm_rank(topo.comm, &rank);
		rank++;
		MPI_Fint fcomm = MPI_Comm_c2f(topo.comm);

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		// TODO: preallocate this
		array<len_t, real_t, 4> yo(fsten.len(0), fsten.len(1), 2, 14);
		int jpn;
		BMG_get_bc(params.per_mask(), &jpn);

		MPI_BMG3_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(),
		                                P.data(), fop.len(0), fop.len(1), fop.len(2),
		                                cop.len(0), cop.len(1), cop.len(2),
		                                nog, ifd, nstencil, BMG_RELAX_SYM,
		                                yo.data(), nog, topo.IGRD(),
		                                ctx->msg_geom.data(), ctx->msg_geom.size(),
		                                ctx->pMSG.data(), ctx->msg_buffer.data(), ctx->msg_buffer.size(),
		                                rank, fcomm, jpn);
	}
}

}}}

#endif
