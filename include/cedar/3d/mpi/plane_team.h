#ifndef CEDAR_3D_MPI_PLANE_TEAM_H
#define CEDAR_3D_MPI_PLANE_TEAM_H

#include <vector>

#include <cedar/2d/mpi/types.h>
#include <cedar/service_manager.h>

namespace cedar { namespace cdr3 { namespace mpi {

class plane_team
{
	using sman_t = service_manager<cdr2::mpi::stypes>;
public:
	std::vector<sman_t*> masters;
	std::vector<std::vector<sman_t*>> workers;
};

}}}


#endif
