#ifndef CEDAR_2D_MPI_INTERP_H
#define CEDAR_2D_MPI_INTERP_H

#include <cedar/2d/mpi/types.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/setup_interp.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_SETUP_interp_OI(int kf, int kc, real_t *so, real_t *ci,
	                                     len_t iif, len_t jjf, len_t iic, len_t jjc,
	                                     int nog, int nogm, len_t *igrd,
	                                     int ifd, int nstncl, int ibc,
	                                     void *halof);
	void BMG_get_bc(int, int*);
}

namespace cedar { namespace cdr2 { namespace mpi {

class interp_f90 : public kernels::interp_add<stypes>
{
	void run(const prolong_op & P,
	         const grid_func & coarse,
	         const grid_func & residual,
	         grid_func & fine) override;

};


class setup_interp_f90 : public kernels::setup_interp<stypes>
{
	void run(const stencil_op<five_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }
	void run(const stencil_op<nine_pt> & fop,
	         const stencil_op<nine_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }

	template<class sten>
	void run_impl(const stencil_op<sten> & fop,
	              const stencil_op<nine_pt> & cop,
	              prolong_op & P)
	{
		int ifd, nstencil, jpn;
		int kf, kc, nog;

		auto & fopd = const_cast<mpi::stencil_op<sten>&>(fop);
		grid_topo & topo = fopd.grid();

		store_fine_op(fopd, P);

		if (stencil_ndirs<sten>::value == 3) {
			ifd = 1;
			nstencil = 3;
		} else {
			ifd = 0;
			nstencil = 5;
		}

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		BMG_get_bc(params->per_mask(), &jpn);

		MPI_BMG2_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), P.data(),
		                                fop.len(0), fop.len(1), cop.len(0), cop.len(1),
		                                nog, nog, topo.IGRD(), ifd, nstencil, jpn, halof);
	}


	inline void store_fine_op(stencil_op<five_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_five = &fop;
		P.fine_is_five = true;
	}

	inline void store_fine_op(stencil_op<nine_pt> & fop,
		                          prolong_op & P)
	{
		P.fine_op_nine = &fop;
		P.fine_is_five = false;
	}
};

}}}

#endif
