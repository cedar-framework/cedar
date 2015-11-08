#ifndef BOXMG_2D_CORE_MPI_GRIDTOPO_H
#define BOXMG_2D_CORE_MPI_GRIDTOPO_H

#include <mpi.h>
#include <array>

#include "boxmg-common.h"

namespace boxmg { namespace bmg2d { namespace mpi {

class grid_topo
{

public:
	using igrd_t = std::shared_ptr<std::vector<len_t>>;
	grid_topo();
	grid_topo(igrd_t grd, int level, int nlevel);

	void grow(int nlevels);
	igrd_t get_igrd() { return igrd; }

	len_t & is(int dim);
	len_t & nlocal(int dim);
	len_t & nglobal(int dim);
	int & nproc(int dim);
	int & coord(int dim);

	const len_t & is (int dim) const;
	const len_t & nlocal(int dim) const;
	const len_t & nglobal(int dim) const;
	const int & nproc(int dim) const;
	const int & coord(int dim) const;

	len_t nproc() const { return nproc_[0]*nproc_[1]; }
	int nlevel() const { return nlvl; }
	int level() const { return lvl; }

	len_t *IGRD() { return igrd->data(); }
	igrd_t igrd_ptr() { return igrd; }

	MPI_Comm comm;
	std::vector<int> dimxfine;
	std::vector<int> dimyfine;

private:
	igrd_t igrd;
	int lvl;
	int nlvl;
	std::array<int, 2> nproc_;
	std::array<int, 2> coord_;
};

using topo_ptr = std::shared_ptr<grid_topo>;

}}}

#endif
