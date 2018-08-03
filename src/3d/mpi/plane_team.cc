#include <cedar/3d/mpi/plane_team.h>

namespace cedar { namespace cdr3 { namespace mpi {

void plane_team::add_worker(sman_t *worker)
{
	if (redist_level >= workers.size())
		workers.push_back({worker});
	else
		workers[redist_level].push_back(worker);

	redist_level = (redist_level + 1) % masters.size();
}

plane_team::sman_t & plane_team::get_master()
{
	return *masters[redist_level];
}


}}}
