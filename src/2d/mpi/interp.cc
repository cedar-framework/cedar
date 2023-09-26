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

using real_t = cedar::real_t;
using len_t = cedar::len_t;

#include <src/2d/ftn/mpi/BMG2_SymStd_interp_add.f90.hpp>

namespace cedar { namespace cdr2 { namespace mpi {

template <typename device>
void interp_f90<device>::run(const prolong_op & P,
                     const grid_func & coarse,
                     const mpi::grid_func & residual,
                     grid_func & fine)
{
		int nstencil, kf, kc, nog;

		auto & Pd = const_cast<prolong_op&>(P);
		grid_func & coarsed = const_cast<grid_func&>(coarse);
		grid_func & res = const_cast<grid_func&>(residual);
		grid_topo & topo = Pd.grid();

                ftl::Buffer<real_t, len_t> fop_data;
		topo_ptr topof;

		if (Pd.fine_is_five) {
			nstencil = 3;
			fop_data = Pd.fine_op_five->to_buffer();
			topof = Pd.fine_op_five->grid_ptr();
		} else {
			nstencil = 5;
			fop_data = Pd.fine_op_nine->to_buffer();
			topof = Pd.fine_op_nine->grid_ptr();
		}

		kc = topo.level() + 1;
		nog = topo.nlevel();
		kf = kc + 1;

                void* halof = services->fortran_handle<mpi::halo_exchange>();

                Pd.template ensure<device>();
                coarsed.template ensure<device>();
                fine.template ensure<device>();

                if (device::is_gpu()) {
                    MPI_BMG2_SymStd_interp_add<cedar::gpu>(
                        kc, kf, nog, fine, coarsed, res,
                        fop_data, nstencil, Pd,
                        coarsed.len(0), coarsed.len(1),
                        fine.len(0), fine.len(1),
                        topof->is(0), topof->is(1), halof);
                } else {
                    MPI_BMG2_SymStd_interp_add(kc, kf, nog,
                                               fine.data(), coarsed.data(), res.data(),
                                               fop_data.data(), nstencil, Pd.data(),
                                               coarsed.len(0), coarsed.len(1),
                                               fine.len(0), fine.len(1),
                                               topof->is(0), topof->is(1),
                                               halof);
                    fine.template mark_dirty<device>();
                }
}

template
void interp_f90<cedar::cpu>::run(const prolong_op&, const grid_func&, const mpi::grid_func&, grid_func&);

template
void interp_f90<cedar::gpu>::run(const prolong_op&, const grid_func&, const mpi::grid_func&, grid_func&);

}}}
