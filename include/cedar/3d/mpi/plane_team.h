#ifndef CEDAR_3D_MPI_PLANE_TEAM_H
#define CEDAR_3D_MPI_PLANE_TEAM_H

#include <vector>
#include <abt.h>

#include <cedar/2d/mpi/types.h>
#include <cedar/service_manager.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_team
{
	using sman_t = service_manager<cdr2::mpi::stypes>;
public:
	plane_team() : redist_level(0) {}

	void add_worker(sman_t *worker);
	sman_t & get_master();

	std::vector<sman_t*> masters;
	std::vector<std::vector<sman_t*>> workers;
	std::vector<ABT_thread> *threads;

protected:
	std::size_t redist_level;
};

}}}


#endif
