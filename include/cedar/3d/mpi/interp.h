#ifndef CEDAR_3D_MPI_INTERP_H
#define CEDAR_3D_MPI_INTERP_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/kernels/interp_add.h>
#include <cedar/kernels/setup_interp.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_interp_OI(int kgf, int kgc, real_t *so, real_t *soc,
	                                     real_t *ci, len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     int nog, int ifd, int nstencil, int irelax, real_t *yo,
	                                     int nogm, len_t *IGRD, int jpn, void *halof);
	void BMG_get_bc(int, int*);
}


namespace cedar { namespace cdr3 { namespace mpi {
class interp_f90 : public kernels::interp_add<stypes>
{
	using prolong_op = cedar::cdr3::mpi::prolong_op;
	using grid_func = cedar::cdr3::mpi::grid_func;
	void run(const prolong_op & P,
	         const grid_func & coarse,
	         const grid_func & residual,
	         grid_func & fine) override;
};


class setup_interp_f90 : public kernels::setup_interp<stypes>
{
	void run(const stencil_op<seven_pt> & fop,
	         const stencil_op<xxvii_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }
	void run(const stencil_op<xxvii_pt> & fop,
	         const stencil_op<xxvii_pt> & cop,
	         prolong_op & P) override { this->run_impl(fop, cop, P); }


	template<class sten>
	void run_impl(const stencil_op<sten> & fop,
	              const stencil_op<xxvii_pt> & cop,
	              prolong_op & P)
	{
		int ifd, nstencil;
		int kc, nog, kf;

		auto & fopd = const_cast<stencil_op<sten>&>(fop);
		auto & copd = const_cast<stencil_op<xxvii_pt>&>(cop);
		grid_topo & topo = fopd.grid();

		store_fine_op(fopd, P);

		nstencil = stencil_ndirs<sten>::value;

		if (std::is_same<sten, seven_pt>::value)
			ifd = 1;
		else
			ifd = 0;

		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		// TODO: preallocate this
		array<real_t, 4> yo(fop.len(0), fop.len(1), 2, 14);
		int jpn;
		BMG_get_bc(params->per_mask(), &jpn);

		MPI_BMG3_SymStd_SETUP_interp_OI(kf, kc, fopd.data(), copd.data(),
		                                P.data(), fop.len(0), fop.len(1), fop.len(2),
		                                cop.len(0), cop.len(1), cop.len(2),
		                                nog, ifd, nstencil, BMG_RELAX_SYM,
		                                yo.data(), nog, topo.IGRD(),
		                                jpn, halof);
	}

	inline void store_fine_op(stencil_op<seven_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_seven = &fop;
		P.fine_is_seven = true;
	}


	inline void store_fine_op(stencil_op<xxvii_pt> & fop,
	                          prolong_op & P)
	{
		P.fine_op_xxvii = &fop;
		P.fine_is_seven = false;
	}
};


}}}

#endif

