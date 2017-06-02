#ifndef CEDAR_2D_CORE_MPI_GRIDFUNC_H
#define CEDAR_2D_CORE_MPI_GRIDFUNC_H


#include <mpi.h>

#include <cedar/2d/grid_func.h>
#include <cedar/mpi/par_object.h>
#include <iostream>

namespace cedar { namespace cdr2 { namespace inter { namespace mpi {
				class prolong_op;
			}
		}
	}
}

namespace cedar { namespace cdr2 { namespace mpi {

class grid_func : public ::cedar::cdr2::grid_func, public par_object
{
public:
	/* template<typename... ArgTypes> */
	/* 	grid_func(ArgTypes&&... args) : ::cedar::cdr2::core::grid_func(std::forward<decltype(args)>(args)...) {} */
	using iadd_t = std::tuple<const cedar::cdr2::inter::mpi::prolong_op&, const grid_func&, const grid_func&>;
	grid_func(){};
	grid_func(topo_ptr grid);
	grid_func(len_t nx, len_t ny);
	using ::cedar::cdr2::grid_func::operator();
	static grid_func like(const grid_func &likeable);
	static grid_func zeros_like(const grid_func &likeable);
	static grid_func ones_like(const grid_func &likeable);
	template<int p> real_t lp_norm() const;
	real_t inf_norm() const;
	grid_func & operator-=(const grid_func & rhs);
	friend grid_func operator-(grid_func lhs, const grid_func &rhs) { return lhs -= rhs; }
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
