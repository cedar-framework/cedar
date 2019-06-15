#include "cedar/2d/ftn/mpi/BMG_workspace_c.h"
#include <cedar/timer.h>
#include <cedar/types.h>
#include <cedar/mpi/grid_topo.h>

#include <cedar/capi.h>

#include <cedar/interface/object.h>

extern "C"
{

	int cedar_topo_create3d(MPI_Comm comm,
	                        unsigned int ngx, unsigned int ngy, unsigned int ngz,
	                        unsigned int nlx[], unsigned int nly[], unsigned int nlz[],
	                        int nprocx, int nprocy, int nprocz,
	                        cedar_topo *newtopo)
	{
		using namespace cedar;

		cedar_object *obj = cedar_object_create(CEDAR_KIND_TOPO);

		int rank;
		MPI_Comm_rank(comm, &rank);

		std::shared_ptr<grid_topo> *grid_ptr;
		auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
		auto grid = std::make_shared<grid_topo>(igrd, 0, 1);

		grid->comm = comm;
		log::set_comm(comm);
		timer_init(comm);

		grid->nproc(0) = nprocx;
		grid->nproc(1) = nprocy;
		grid->nproc(2) = nprocz;

		grid->coord(0) = (rank % (nprocx * nprocy)) % nprocx;
		grid->coord(1) = (rank % (nprocx * nprocy)) / nprocx;
		grid->coord(2) = rank / (nprocx * nprocy);

		grid->is(0) = 1; grid->is(1) = 1; grid->is(2) = 1;
		for (auto i : range(grid->coord(0))) {
			grid->is(0) += nlx[i];
		}

		for (auto j : range(grid->coord(1))) {
			grid->is(1) += nly[j];
		}

		for (auto k : range(grid->coord(2))) {
			grid->is(2) += nlz[k];
		}

		grid->nglobal(0) = ngx + 2;
		grid->nglobal(1) = ngy + 2;
		grid->nglobal(2) = ngz + 2;

		grid->nlocal(0) = nlx[grid->coord(0)] + 2;
		grid->nlocal(1) = nly[grid->coord(1)] + 2;
		grid->nlocal(2) = nlz[grid->coord(2)] + 2;

		grid->dimxfine.resize(nprocx);
		grid->dimyfine.resize(nprocy);
		grid->dimzfine.resize(nprocz);
		for (auto i: range(nprocx)) {
			grid->dimxfine[i] = nlx[i];
		}

		for (auto j: range(nprocy)) {
			grid->dimyfine[j] = nly[j];
		}

		for (auto k: range(nprocz)) {
			grid->dimzfine[k] = nlz[k];
		}

		grid_ptr = new std::shared_ptr<grid_topo>(std::move(grid));
		obj->ptr = reinterpret_cast<void*>(grid_ptr);
		*newtopo = obj->handle;

		return CEDAR_SUCCESS;
	}
}
