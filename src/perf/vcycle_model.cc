#include <cmath>

#include <boxmg/perf/vcycle_model.h>


using namespace boxmg;

vcycle_model::vcycle_model(short nd, int v1, int v2) : nd(nd), v1(v1), v2(v2)
{
	if (nd == 2) {
		nc = 4;
		ns = 9;
	} else {
		nc = 8;
		ns = 27;
	}
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


float vcycle_model::tsmooth(int lvl)
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


float vcycle_model::tresidual(int lvl)
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


float vcycle_model::trestrict(int lvl)
{
	float time = 0;
	float prod = 1;
	for (int i = 0; i < nd; i++) {
		prod *= grid(lvl).nlocal(i);
	}

	time += 2*ns*prod*tc;

	return time;
}


float vcycle_model::tinterp(int lvl)
{
	float time = 0;

	if (nd == 2) {
		time += grid(lvl).nlocal(0)*grid(lvl).nlocal(1) * tc;
		time += 20 * grid(lvl-1).nlocal(0)*grid(lvl-1).nlocal(1) * tc;
		time += 6 * (grid(lvl-1).nlocal(0) + grid(lvl-1).nlocal(1)) * tc;

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
		time += (prodf + 60 * prodc + 15 * grid(lvl-1).nlocal(0) * grid(lvl-1).nlocal(3)) * tc;
		time += (6*grid(lvl-1).nlocal(1)*grid(lvl-1).nlocal(2) + grid(lvl-1).nlocal(2)) * tc;
		time += 6*ts + 2*sigmaf * tw;
	}

	return time;
}


float vcycle_model::tcgsolve()
{
	float time = 0;
	int p = grid(0).nproc();

	time += cg_perf->time();
	time += std::ceil(std::log2(p))*ts;
	float prod = 1;
	for (int i = 0; i < nd; i++) {
		prod *= grid(0).nglobal(i);
	}
	time += prod*(1 + std::ceil(std::log2(p)))*tw;

	return time;
}


float vcycle_model::time()
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
