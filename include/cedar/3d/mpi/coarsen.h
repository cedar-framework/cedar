#ifndef CEDAR_3D_MPI_COARSEN_H
#define CEDAR_3D_MPI_COARSEN_H

#include <cedar/2d/ftn/mpi/BMG_parameters_c.h>
#include <cedar/kernels/coarsen_op.h>
#include <cedar/3d/mpi/types.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG3_SymStd_SETUP_ITLI27_ex(int kgf, int kgc, real_t *so, real_t *soc, real_t *ci,
	                                     len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     len_t iGs, len_t jGs, len_t kGs,
	                                     int nogm, len_t *IGRD,
	                                     int nproc, void *halof);
	void MPI_BMG3_SymStd_SETUP_ITLI07_ex(int kgf, int kgc, real_t *so, real_t *soc, real_t *ci,
	                                     len_t iif, len_t jjf, len_t kkf,
	                                     len_t iic, len_t jjc, len_t kkc,
	                                     len_t iGs, len_t jGs, len_t kGs,
	                                     void *halof);
}

namespace cedar { namespace cdr3 { namespace mpi {

class galerkin : public kernels::coarsen_op<stypes>
{
	void run(const prolong_op & P,
	         const stencil_op<seven_pt> & fop,
	         stencil_op<xxvii_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}
	void run(const prolong_op & P,
	         const stencil_op<xxvii_pt> & fop,
	         stencil_op<xxvii_pt> & cop) override
	{
		this->run_impl(P, fop, cop);
	}


	template<class sten>
	void run_impl(const prolong_op & P,
	              const stencil_op<sten> & fop,
	              stencil_op<xxvii_pt> & cop)
	{
		int kf, kc, nog;
		prolong_op & Pd = const_cast<prolong_op&>(P);
		auto & fopd = const_cast<stencil_op<sten>&>(fop);
		grid_topo & topo = fopd.grid();
		kc = topo.level();
		nog = topo.nlevel();
		kf = kc + 1;

		if (std::is_same<sten, seven_pt>::value) {
			MPI_BMG3_SymStd_SETUP_ITLI07_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
			                                fop.len(0), fop.len(1), fop.len(2),
			                                cop.len(0), cop.len(1), cop.len(2),
			                                topo.is(0), topo.is(1), topo.is(2),
			                                halof);
		} else {
			MPI_BMG3_SymStd_SETUP_ITLI27_ex(kf, kc, fopd.data(), cop.data(), Pd.data(),
			                                fop.len(0), fop.len(1), fop.len(2),
			                                cop.len(0), cop.len(1), cop.len(2),
			                                topo.is(0), topo.is(1), topo.is(2),
			                                nog, topo.IGRD(), topo.nproc(), halof);
		}
	}
};

}}}

#endif
