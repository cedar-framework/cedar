#include <cedar/2d/ftn/BMG_parameters_c.h>
#include <cedar/2d/mpi/interp.h>

extern "C" {
	using namespace cedar;
	void MPI_BMG2_SymStd_interp_add(len_t kc, len_t kf, len_t nog,
	                                real_t *Q, real_t *QC, real_t *RES, real_t *SO,
	                                int nstncl, real_t *CI, len_t iic, len_t jjc,
	                                len_t iif, len_t jjf, len_t iGs, len_t jGs,
	                                void *halof);
}

namespace cedar { namespace cdr2 { namespace mpi {

void interp_f90::run(const prolong_op & P,
                     const grid_func & coarse,
                     const mpi::grid_func & residual,
                     grid_func & fine)
{
		int nstencil, kf, kc, nog;

		auto & Pd = const_cast<prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);
		grid_topo & topo = Pd.grid();

		real_t * fop_data;
		topo_ptr topof;

		if (Pd.fine_is_five) {
			nstencil = 3;
			fop_data = Pd.fine_op_five->data();
			topof = Pd.fine_op_five->grid_ptr();
		} else {
			nstencil = 5;
			fop_data = Pd.fine_op_nine->data();
			topof = Pd.fine_op_nine->grid_ptr();
		}

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

		MPI_BMG2_SymStd_interp_add(kc, kf, nog,
		                           fine.data(), coarsed.data(), res.data(),
		                           fop_data, nstencil,
		                           Pd.data(),
		                           coarsed.len(0), coarsed.len(1),
		                           fine.len(0), fine.len(1),
		                           topof->is(0), topof->is(1),
		                           services->template fortran_handle<halo_exchange>());
}

}}}
