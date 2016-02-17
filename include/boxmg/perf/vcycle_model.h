#ifndef BOXMG_VCYCLE_MODEL_H
#define BOXMG_VCYCLE_MODEL_H

#include <vector>
#include <memory>

#include <boxmg/mpi/grid_topo.h>
#include <boxmg/perf/perf_model.h>

namespace boxmg {

class vcycle_model : public perf_model
{
public:
	vcycle_model(short nd, int v1=2, int v2=1);
	void add_level(topo_ptr grid);
	grid_topo & grid(int lvl);
	int ngrids() { return grids.size(); }
	float tsmooth(int lvl);
	float tresidual(int lvl);
	float trestrict(int lvl);
	float tinterp(int lvl);
	float tcgsolve();
	virtual float time();

protected:
	short nd;
	std::vector<topo_ptr> grids;
	int nc, ns;
	int v1, v2;
	std::shared_ptr<perf_model> cg_perf;
};

}

#endif
