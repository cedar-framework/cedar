#ifndef BOXMG_2D_KERNEL_MPI_REGISTRY_H
#define BOXMG_2D_KERNEL_MPI_REGISTRY_H


#include "boxmg/types.h"
#include "boxmg/kernel_registry.h"
#include "boxmg/cycle/types.h"
#include "boxmg/kernel_name.h"

#include "boxmg/2d/relax_stencil.h"
#include "boxmg/2d/mpi/grid_func.h"
#include "boxmg/2d/mpi/grid_topo.h"
#include "boxmg/2d/mpi/grid_func.h"


namespace boxmg { namespace bmg2d {
			class solver;
}}

namespace boxmg { namespace bmg2d { namespace mpi {
			class stencil_op;
		}
		namespace inter {
			namespace mpi {
				class prolong_op;
			}
		}
	}
}


namespace boxmg { namespace bmg2d { namespace kernel { namespace mpi {
namespace mpi = boxmg::bmg2d::mpi;
class registry : public kernel_registry<mpi::stencil_op, relax_stencil, inter::mpi::prolong_op, mpi::grid_func>
{
public:
	void setup_nog(mpi::grid_topo &topo,
	               len_t min_coarse, int *nog)
	{
		active.run(kernel_name::setup_nog, topo,
		           static_cast<len_t>(min_coarse),
		           static_cast<int*>(nog));
	}


	void halo_setup(mpi::grid_topo &topo,
	                void **halo_ctx)
	{
		active.run(kernel_name::halo_setup,
		           static_cast<mpi::grid_topo&>(topo),
		           static_cast<void**>(halo_ctx));
	}


	void halo_exchange(mpi::grid_func &f)
	{
		active.run(kernel_name::halo_exchange,
		           static_cast<mpi::grid_func&>(f));
	}


	void halo_exchange(const mpi::grid_func &f, void *halo_ctx)
	{
		auto & fd = const_cast<mpi::grid_func&>(f);
		fd.halo_ctx = halo_ctx;
		active.run(kernel_name::halo_exchange,
		           static_cast<mpi::grid_func&>(fd));

	}


	void halo_stencil_exchange(mpi::stencil_op & so)
	{
		active.run(kernel_name::halo_stencil_exchange,
		           so);
	}


	void setup_cg_boxmg(const mpi::stencil_op & so,
	                    std::shared_ptr<solver> *bmg)
	{
		active.run(kernel_name::setup_cg_boxmg, so,
		           static_cast<std::shared_ptr<solver>*>(bmg));

	}


	void solve_cg_boxmg(const solver &bmg,
	                    grid_func &x,
	                    const grid_func &b)
	{
		active.run(kernel_name::solve_cg_boxmg, bmg, x, b);
	}


	void matvec(const mpi::stencil_op & so,
	            const grid_func &x,
	            grid_func &b)
	{
		active.run(kernel_name::matvec, so, x, b);
	}
};

}}}}


#endif
