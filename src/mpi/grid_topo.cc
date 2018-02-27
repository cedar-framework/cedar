#include <cedar/2d/ftn/mpi/BMG_workspace_c.h>

#include <cedar/mpi/grid_topo.h>

using namespace cedar;

grid_topo::grid_topo() : igrd(nullptr) {}

grid_topo::grid_topo(igrd_t grd, int level, int nlevel) : igrd(grd), lvl(level), nlvl(nlevel) {}

len_t & grid_topo::is(int dim)
{
	if (dim == 0) {
		return (*igrd)[(idL_BMG_Icoord-1)*nlvl + lvl];
	} else if (dim == 1) {
		return (*igrd)[(idL_BMG_Jcoord-1)*nlvl + lvl];
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_Kcoord-1)*nlvl + lvl];
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
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_Kcoord-1)*nlvl + lvl];
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
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_NLz-1)*nlvl + lvl];
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
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_NLz-1)*nlvl + lvl];
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
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_NGz-1)*nlvl + lvl];
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
	} else if (dim == 2) {
		return (*igrd)[(idL_BMG_NGz-1)*nlvl + lvl];
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


namespace cedar {
	std::ostream & operator<<(std::ostream & os, const grid_topo & obj)
	{
		os << "===== grid_topo =====\n";

		os << "nproc:   " << obj.nproc(0) << " " << obj.nproc(1) << " " << obj.nproc(2) << '\n';
		os << "coord:   " << obj.coord(0) << " " << obj.coord(1) << " " << obj.coord(2) << '\n';
		os << "nlocal:  " << obj.nlocal(0) << " " << obj.nlocal(1) << " " << obj.nlocal(2) << '\n';
		os << "nglobal: " << obj.nglobal(0) << " " << obj.nglobal(1) << " " << obj.nglobal(2) << '\n';
		os << "is:      " << obj.is(0) << " " << obj.is(1) << " " << obj.is(2) << '\n';

		return os;
	}
}
