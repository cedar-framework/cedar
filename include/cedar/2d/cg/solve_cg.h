#ifndef CEDAR_2D_KERNEL_SOLVE_CG_H
#define CEDAR_2D_KERNEL_SOLVE_CG_H

#include <cedar/kernel_params.h>
#include <cedar/2d/grid_func.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/halo_exchanger_base.h>
#include <cedar/2d/mpi/msg_exchanger.h>

namespace cedar { namespace cdr2 {
		template<class> class solver;
		enum class nine_pt;
}}

namespace cedar { namespace cdr2 { namespace kernel {

namespace impls
{
	void fortran_solve_cg(const kernel_params & params, grid_func & x,
	                      const grid_func & b,
	                      const grid_func & ABD,
	                      real_t *bbd);


	void mpi_solve_cg_lu(const kernel_params & params,
	                     mpi::msg_exchanger *halof,
	                     mpi::grid_func &x,
	                     const mpi::grid_func &b,
	                     const mpi::grid_func & ABD,
	                     real_t *bbd);


	void solve_cg_boxmg(const kernel_params & params,
	                    mpi::msg_exchanger *halof,
	                    const solver<nine_pt> & cg_solver,
	                    mpi::grid_func &x,
	                    const mpi::grid_func &b);
}

}}}

#endif
