#ifndef BOXMG_2D_CORE_MPI_GRIDFUNC_H
#define BOXMG_2D_CORE_MPI_GRIDFUNC_H


#include <mpi.h>

#include <boxmg/2d/grid_func.h>
#include <boxmg/mpi/par_object.h>
#include <iostream>

namespace boxmg { namespace bmg2d { namespace inter { namespace mpi {
				class prolong_op;
			}
		}
	}
}

namespace boxmg { namespace bmg2d { namespace mpi {

class grid_func : public ::boxmg::bmg2d::grid_func, public par_object
{
public:
	/* template<typename... ArgTypes> */
	/* 	grid_func(ArgTypes&&... args) : ::boxmg::bmg2d::core::grid_func(std::forward<decltype(args)>(args)...) {} */
	using iadd_t = std::tuple<const boxmg::bmg2d::inter::mpi::prolong_op&, const grid_func&, const grid_func&>;
	grid_func(){};
	grid_func(topo_ptr grid);
	grid_func(len_t nx, len_t ny);
	using ::boxmg::bmg2d::grid_func::operator();
	static grid_func like(const grid_func &likeable);
	static grid_func zeros_like(const grid_func &likeable);
	static grid_func ones_like(const grid_func &likeable);
	template<int p> real_t lp_norm() const;
	real_t inf_norm() const;
	grid_func & operator-=(const grid_func & rhs);
	friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
	grid_func & operator+=(iadd_t interp_add_package);
	friend std::ostream& operator<< (std::ostream& os, const grid_func &obj);
};


template<int p> real_t grid_func::lp_norm() const
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

}}}

#endif
