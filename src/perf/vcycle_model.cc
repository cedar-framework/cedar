#include <cmath>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <cedar/perf/vcycle_model.h>


using namespace cedar;

vcycle_model::vcycle_model(short nd, int v1, int v2) : nd(nd), v1(v1), v2(v2), isleaf(true)
{
	if (nd == 2) {
		nc = 4;
		ns = 9;
	} else {
		nc = 8;
		ns = 27;
	}

	for (auto i : range(nd)) {
		nb[i] = 1;
	}
}


int & vcycle_model::nblocks(int dim)
{
	return this->nb[dim];
}


int vcycle_model::nblocks(int dim) const
{
	return this->nb[dim];
}


void vcycle_model::add_level(topo_ptr grid)
{
	grids.push_back(grid);
}


grid_topo & vcycle_model::grid(int i)
{
	int ind;
	if (i < 0) ind = (i + 1)*-1;
	else ind = grids.size() - i - 1;
	#ifdef DEBUG
	return *grids.at(ind);
	#else
	return *grids[ind];
	#endif
}


const grid_topo & vcycle_model::grid(int i) const
{
	int ind;
	if (i < 0) ind = (i + 1)*-1;
	else ind = grids.size() - i - 1;
	#ifdef DEBUG
	return *grids.at(ind);
	#else
	return *grids[ind];
	#endif
}


topo_ptr vcycle_model::grid_ptr(int i)
{
	int ind;
	if (i < 0) ind = (i + 1)*-1;
	else ind = grids.size() - i - 1;
	#ifdef DEBUG
	return grids.at(ind);
	#else
	return grids[ind];
	#endif
}


float vcycle_model::tsmooth(int lvl) const
{
	float time = 0;
	float prod = 1;
	float sigma = 0;
	for (int i = 0; i < nd; i++) {
		prod *= grid(lvl).nlocal(i);
		sigma += grid(lvl).nlocal(i);
	}
	time += 2 * ns * prod * (v1+v2) * tc;
	time += (v1 + v2)*(2*nd*ts + 2*sigma*tw)*nc;

	return time;
}


float vcycle_model::tresidual(int lvl) const
{
	float time = 0;
	float prod = 1;
	float sigma = 0;
	for (int i = 0; i < nd; i++) {
		prod *= grid(lvl).nlocal(i);
		sigma += grid(lvl).nlocal(i);
	}

	time += 2*ns*prod*tc;
	time += 2*nd*ts + 2*sigma*tw;

	return time;
}


float vcycle_model::trestrict(int lvl) const
{
	float time = 0;
	float prod = 1;
	for (int i = 0; i < nd; i++) {
		prod *= grid(lvl).nlocal(i);
	}

	time += 2*ns*prod*tc;

	return time;
}


float vcycle_model::tinterp(int lvl) const
{
	float time = 0;
	float tc_interp =  2.7828045360610174*tc;

	if (nd == 2) {
		time += grid(lvl).nlocal(0)*grid(lvl).nlocal(1) * tc_interp;
		time += 20 * grid(lvl-1).nlocal(0)*grid(lvl-1).nlocal(1) * tc_interp;
		time += 6 * (grid(lvl-1).nlocal(0) + grid(lvl-1).nlocal(1)) * tc_interp;

		time += 4 * ts;
		time += 2 * (grid(lvl).nlocal(0) + grid(lvl).nlocal(1)) * tw;
	} else {
		float prodf = 1;
		float prodc = 1;
		float sigmaf = 0;
		for (int i = 0; i < nd; i++) {
			prodf *= grid(lvl).nlocal(i);
			prodc *= grid(lvl-1).nlocal(i);
			sigmaf += grid(lvl).nlocal(i);
		}
		time += (prodf + 60 * prodc + 15 * grid(lvl-1).nlocal(0) * grid(lvl-1).nlocal(3)) * tc_interp;
		time += (6*grid(lvl-1).nlocal(1)*grid(lvl-1).nlocal(2) + grid(lvl-1).nlocal(2)) * tc_interp;
		time += 6*ts + 2*sigmaf * tw;
	}

	return time;
}


float vcycle_model::tcgsolve() const
{
	float time = 0;

	if (!isleaf) {
		// int p = grid(0).nproc();
		int gather_size = 1;
		len_t cg_size = 1;
		for (auto i : range(nd)) {
			gather_size *= std::ceil(static_cast<double>(grid(0).nproc(i)) / static_cast<double>(nb[i]));
			cg_size *= std::ceil(static_cast<double>(grid(0).nglobal(i))/static_cast<double>(nb[i]));
		}

		// allgather
		time += cg_perf->time();
		time += std::ceil(std::log2(gather_size))*ts;
		time += cg_size*(1 + std::ceil(std::log2(gather_size)))*tw;
		// time += cg_size*(gather_size-1)/gather_size*tw;

		// gather scatter
		// time += 2*std::ceil(std::log2(gather_size))*ts;
		// time += cg_size*tw;

		// time += 2*std::ceil(std::log2(gather_size))*ts;
		// time += cg_size*tw;
	}

	return time;
}



float vcycle_model::agglom() const
{
	float time = 0;

	if (!isleaf) {
		// int p = grid(0).nproc();
		int gather_size = 1;
		len_t cg_size = 1;
		for (auto i : range(nd)) {
			gather_size *= std::ceil(static_cast<double>(grid(0).nproc(i)) / static_cast<double>(nb[i]));
			cg_size *= std::ceil(static_cast<double>(grid(0).nglobal(i))/static_cast<double>(nb[i]));
		}

		// allgather
		time += std::ceil(std::log2(gather_size))*ts;
		time += cg_size*(1 + std::ceil(std::log2(gather_size)))*tw;
		// time += cg_size*(gather_size-1)/gather_size*tw;

		// gather scatter
		// time += 2*std::ceil(std::log2(gather_size))*ts;
		// time += cg_size*tw;

		// time += 2*std::ceil(std::log2(gather_size))*ts;
		// time += cg_size*tw;
	}

	return time;
}


float vcycle_model::time() const
{
	float time = 0;

	time += tcgsolve();
	for (int i = 1; i < ngrids(); i++) {
		time += tsmooth(i);
		time += tresidual(i);
		time += trestrict(i);
		time += tinterp(i);
	}

	return time;
}


void vcycle_model::set_cgperf(std::shared_ptr<perf_model> mod)
{
	cg_perf = mod;
	isleaf = false;
}


std::shared_ptr<perf_model> vcycle_model::get_cgperf()
{
	return cg_perf;
}


int vcycle_model::nproc() const
{
	return grid(0).nproc();
}


void vcycle_model::save_levels()
{
	std::ofstream smooth_file;
	std::ofstream residual_file;
	std::ofstream interp_file;
	std::ofstream restrict_file;
	smooth_file.open("smooth.txt", std::ios::out | std::ios::trunc | std::ios::binary);
	residual_file.open("residual.txt", std::ios::out | std::ios::trunc | std::ios::binary);
	interp_file.open("interp.txt", std::ios::out | std::ios::trunc | std::ios::binary);
	restrict_file.open("restrict.txt", std::ios::out | std::ios::trunc | std::ios::binary);
	for (auto i = 1; i < ngrids(); i++) {
		smooth_file << i << " " << tsmooth(i) << '\n';
		residual_file << i << " " << tresidual(i) << '\n';
		interp_file << i << " " << tinterp(i) << '\n';
		restrict_file << i << " " << trestrict(i) << '\n';
	}
	smooth_file.close();
	residual_file.close();
	interp_file.close();
	restrict_file.close();
}


void vcycle_model::rep(std::ostream & os) const
{
	os << "======== vcycle model ========" << '\n';
	os << "nproc:       " << grid(-1).nproc(0) << " " << grid(-1).nproc(1) << '\n';
	os << "nblock:      " << nblocks(0) << " " << nblocks(1) << '\n';
	os << "local size:  " << grid(-1).nlocal(0) << " x " << grid(-1).nlocal(1) << '\n';
	os << "global size: " << grid(-1).nglobal(0) << " x " << grid(-1).nglobal(1) << '\n';
	os << "coarse size: " << grid(0).nglobal(0) << " x " << grid(0).nglobal(1) << '\n';
	os << "nlevel:      " << ngrids() << '\n';
	os << "cgtime:      " << tcgsolve() << '\n';
	os << "time:        " << time() << '\n';
	os << '\n';
	if (!isleaf) {
		os << *cg_perf;
	}
}


void vcycle_model::recur_times(boost::property_tree::ptree & children) const
{
	using namespace boost::property_tree;

	ptree child;

	float smooth_time = 0;
	float residual_time = 0;
	float interp_time = 0;
	float restrict_time = 0;
	for (int i=1; i < ngrids(); i++) {
		smooth_time += tsmooth(i);
		residual_time += tresidual(i);
		interp_time += tinterp(i);
		restrict_time += trestrict(i);
	}

	child.put("relaxation", smooth_time*10);
	child.put("residual", residual_time*10);
	child.put("interp-add", interp_time*10);
	child.put("solve", time()*10);
	child.put("restrict", restrict_time*10);
	child.put("agglomerate", agglom()*10);

	if (!isleaf) {
		(*cg_perf).recur_times(children);
	}

	children.push_back(std::make_pair("",child));
}


void vcycle_model::save_times(std::string iname, std::string oname) const
{
	using namespace boost::property_tree;

	ptree pt, children;

	json_parser::read_json(iname, pt);

	recur_times(children);

	pt.add_child("model", children);
	json_parser::write_json(oname, pt);
}
