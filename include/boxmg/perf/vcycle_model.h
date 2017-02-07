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
	const grid_topo & grid(int lvl) const;
	topo_ptr grid_ptr(int lvl);
	int ngrids() const { return grids.size(); }
	float tsmooth(int lvl) const;
	float tresidual(int lvl) const;
	float trestrict(int lvl) const;
	float tinterp(int lvl) const;
	float tcgsolve() const;
	float agglom() const;
	void set_cgperf(std::shared_ptr<perf_model> mod);
	std::shared_ptr<perf_model> get_cgperf();
	virtual float time() const;
	void save_levels();
	int & nblocks(int dim);
	int nblocks(int dim) const;
	virtual void rep(std::ostream & os) const;
	virtual int nproc() const;
	virtual void recur_times(boost::property_tree::ptree&) const;
	void save_times(std::string iname, std::string oname) const;

protected:
	short nd;
	std::vector<topo_ptr> grids;
	int nc, ns;
	int v1, v2;
	std::shared_ptr<perf_model> cg_perf;
	std::array<int,3> nb;
	bool isleaf;
};

}

#endif
