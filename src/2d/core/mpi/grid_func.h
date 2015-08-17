#ifndef BOXMG_2D_CORE_MPI_GRIDFUNC_H
#define BOXMG_2D_CORE_MPI_GRIDFUNC_H


#include <mpi.h>

#include "core/grid_func.h"
#include "grid_topo.h"

namespace boxmg { namespace bmg2d { namespace core { namespace mpi {

class GridFunc : public ::boxmg::bmg2d::core::GridFunc
{
public:
	/* template<typename... ArgTypes> */
	/* 	GridFunc(ArgTypes&&... args) : ::boxmg::bmg2d::core::GridFunc(std::forward<decltype(args)>(args)...) {} */
	GridFunc(topo_ptr grid);

	MPI_Comm comm;
	GridTopo & grid() { return *grid_; }
	const GridTopo & grid() const { return *grid_;}
	topo_ptr grid_ptr() const { return grid_; }
	void *halo_ctx;
	static GridFunc zeros_like(const GridFunc &likeable);
	static GridFunc ones_like(const GridFunc &likeable);
	template<int p> real_t lp_norm() const;
	real_t inf_norm() const;
	GridFunc & operator-=(const GridFunc & rhs);
	friend GridFunc operator-(GridFunc lhs, const GridFunc &rhs) { return lhs -= rhs; }

protected:
	topo_ptr grid_;
};


template<int p> real_t GridFunc::lp_norm() const
{
	real_t result = 0;

	for (auto j : this->range(1)) {
		for (auto i : this->range(0)) {
			result += std::pow((*this)(i,j), p);
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, &result, 1, MPI_DOUBLE, MPI_SUM, grid_->comm);

	return std::pow(result, 1./p);
}

}}}}

#endif
