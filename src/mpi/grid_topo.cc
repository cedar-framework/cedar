#include <boxmg/2d/ftn/mpi/BMG_workspace_c.h>

#include <boxmg/mpi/grid_topo.h>

using namespace boxmg;

grid_topo::grid_topo() : igrd(nullptr) {}

grid_topo::grid_topo(igrd_t grd, int level, int nlevel) : igrd(grd), lvl(level), nlvl(nlevel) {}

len_t & grid_topo::is(int dim)
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


const len_t & grid_topo::is(int dim) const
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


len_t & grid_topo::nlocal(int dim)
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


const len_t & grid_topo::nlocal(int dim) const
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


len_t & grid_topo::nglobal(int dim)
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


const len_t & grid_topo::nglobal(int dim) const
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


int & grid_topo::nproc(int dim)
{
	return nproc_[dim];
}


const int & grid_topo::nproc(int dim) const
{
	return nproc_[dim];
}


int & grid_topo::coord(int dim)
{
	return coord_[dim];
}


const int & grid_topo::coord(int dim) const
{
	return coord_[dim];
}


void grid_topo::grow(int nlevels)
{
	std::vector<len_t> tmp = (*igrd);

	igrd->resize(nlevels*NBMG_pIGRD);

	for (auto i: range(NBMG_pIGRD)) {
		(*igrd)[i*nlevels + nlevels-1] = tmp[i*nlvl + lvl];
	}

	nlvl = nlevels;
	lvl = nlevels - 1;
}
