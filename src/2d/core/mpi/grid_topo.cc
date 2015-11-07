#include "kernel/fortran/mpi/BMG_workspace_c.h"

#include "grid_topo.h"

using namespace boxmg;
using namespace boxmg::bmg2d::mpi;

GridTopo::GridTopo() : igrd(nullptr) {}

GridTopo::GridTopo(igrd_t grd, int level, int nlevel) : igrd(grd), lvl(level), nlvl(nlevel) {}

len_t & GridTopo::is(int dim)
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_Icoord-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_Jcoord-1)*nlvl + lvl];
	} else {
		log::error << "Invalid dimension" << std::endl;
		return (*igrd)[(idL_BMG_Kcoord-1)*nlvl + lvl];
	}
}


const len_t & GridTopo::is(int dim) const
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_Icoord-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_Jcoord-1)*nlvl + lvl];
	} else {
		log::error << "Invalid dimension" << std::endl;
		return (*igrd)[(idL_BMG_Kcoord-1)*nlvl + lvl];
	}
}


len_t & GridTopo::nlocal(int dim)
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_NLx-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_NLy-1)*nlvl + lvl];
	} else {
		log::error << "invalid dimension" << std::endl;
		return (*igrd)[(idL_BMG_NLz-1)*nlvl + lvl];
	}
}


const len_t & GridTopo::nlocal(int dim) const
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_NLx-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_NLy-1)*nlvl + lvl];
	} else {
		log::error << "invalid dimension" << std::endl;
		return (*igrd)[(idL_BMG_NLz-1)*nlvl + lvl];
	}
}


len_t & GridTopo::nglobal(int dim)
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_NGx-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_NGy-1)*nlvl + lvl];
	} else {
		log::error << "invalid dimension" << std::endl;
		return (*igrd)[0];
	}
}


const len_t & GridTopo::nglobal(int dim) const
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_NGx-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_NGy-1)*nlvl + lvl];
	} else {
		log::error << "invalid dimension" << std::endl;
		return (*igrd)[0];
	}
}


int & GridTopo::nproc(int dim)
{
	return nproc_[dim];
}


const int & GridTopo::nproc(int dim) const
{
	return nproc_[dim];
}


int & GridTopo::coord(int dim)
{
	return coord_[dim];
}


const int & GridTopo::coord(int dim) const
{
	return coord_[dim];
}


void GridTopo::grow(int nlevels)
{
	std::vector<len_t> tmp = (*igrd);

	igrd->resize(nlevels*NBMG_pIGRD);

	for (auto i: range(NBMG_pIGRD)) {
		(*igrd)[i*nlevels + nlevels-1] = tmp[i*nlvl + lvl];
	}

	nlvl = nlevels;
	lvl = nlevels - 1;
}
