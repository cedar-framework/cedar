#ifndef CEDAR_2D_MPI_REDIST_SOLVER_H
#define CEDAR_2D_MPI_REDIST_SOLVER_H

#include <mpi.h>
#include <array>
#include <cassert>

#include <cedar/config.h>
#include <cedar/mpi/redist_comms.h>
#include <cedar/2d/mpi/stencil_op.h>
#include <cedar/2d/mpi/msg_exchanger.h>
#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/array.h>
#include <cedar/mpi/block_partition.h>
#include <cedar/util/timer.h>
#include <cedar/2d/mpi/kernel_manager.h>


extern "C" {
	void MSG_pause(MPI_Fint *msg_comm);
	void MSG_play(MPI_Fint msg_comm);
}


namespace cedar { namespace cdr2 { namespace mpi {


/**
 Coarse-grid redistribution solver.

 Provides coarse-grid redistribution through a generic solver
 interface.  This works by wrapping a solver object distributed on
 fewer cores than the callers current data distribution.

 @note This class does not currently support five_pt operators.
 */
template<class inner_solver>
class redist_solver
{
public:
	/**
	   Performs redistribution of the grid topology and stencil
	   operator given the destination 2D distribution.  This also
	   runs the setup phase of the redistributed solver.

	   @param[in] so The stencil operator to redistribute
	   @param[in] conf The config object to use for this solver
	   @param[in] nblock The destination 2D distribution
	*/
	redist_solver(const stencil_op<nine_pt> & so,
	              service_manager<stypes> * services,
	              std::shared_ptr<config> conf,
	              std::array<int, 2> nblock);
	/**
	   Runs the redistributed solve phase
	   @param[in] b rhs
	   @param[out] x solution
	*/
	void solve(grid_func & x, const grid_func & b);

	inner_solver & get_inner() { return *slv; }

	bool isactive() {
		return (redundant and active) or (not redundant and block_id == 0);
	}

protected:
	bool redundant; /** Flag for whether redistribution is performed redundantly */
	std::unique_ptr<inner_solver> slv; /** Redistributed solver object. */
	redist_comms rcomms; /** Useful communicators used in redistribution */
	int block_id; /** id within a block */
	int block_num; /** which block */
	std::array<int,2> nblock; /** number of processors in each dimension for redistributed solve */
	bool active; /** Flag for whether the current processor is active in the redistribution */
	int recv_id; /** Where the current processor will be receiving data in the fixup phase */
	std::vector<int> send_ids; /** Where the current processor will send data in the fixup phase */
	array<len_t, 1> nbx; /** number of d.o.f. for each processor in my block */
	array<len_t, 1> nby; /** number of d.o.f. for each processor in my block */
	MPI_Fint msg_comm;
	typename inner_solver::grid_func b_redist; /** The redistributed rhs */
	typename inner_solver::grid_func x_redist; /** The redistributed solution */
	std::unique_ptr<typename inner_solver::stencil_op> so_redist; /** The redistributed operator */
	service_manager<stypes> * services;

	/** Redistributes the processor grid.

	    @param[in] fine_topo The source topology to redistribute
	    @param[in] dimx number of grid points on each processor. shape(proc_coord_x, level)
	    @param[in] dimy number of grid points on each processor. shape(proc_coord_y, level)
	    @returns The redistrubted processor grid topology
	*/
	std::shared_ptr<grid_topo> redist_topo(const grid_topo & fine_topo,
	                                       aarray<int, len_t, 2> & dimx,
	                                       aarray<int, len_t, 2> & dimy);
	/** Redistributes the operator.

	    @param[in] so Source operator to redistribute
	    @param[in] topo Redistributed processor grid topology
	*/
	void redist_operator(const stencil_op<nine_pt> & so, topo_ptr topo);
	/** Gathers the rhs

	    @param[in] b The rhs to redistribute
	*/
	void gather_rhs(const grid_func & b);
	/** Scatters the redistributed solution

	    @param[out] x Grid function to receieve the redistributed solution
	*/
	void scatter_sol(grid_func & x);
};


template <class T>
	std::unique_ptr<T> create_operator(topo_ptr topo);

}}}

#endif
