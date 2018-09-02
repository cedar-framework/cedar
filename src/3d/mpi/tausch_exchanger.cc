#include <cedar/3d/mpi/tausch_exchanger.h>


extern "C" {
	using namespace cedar;
	void BMG3_SymStd_SETUP_Tausch(int nog, len_t *dimx, len_t *dimy, len_t *dimz,
	                              len_t *dimxfine, len_t *dimyfine, len_t *dimzfine,
	                              int nproci, int nprocj, int nprock);
}


namespace cedar { namespace cdr3 { namespace mpi {
line_pkg::line_pkg(grid_topo & topo):
	datadist{{array<len_t, 2>(2, topo.nproc(0)),
			array<len_t, 2>(2, topo.nproc(1)),
			array<len_t,2>(2, topo.nproc(2))}}
{
}


void tausch_exchanger::setup(std::vector<topo_ptr> topos)
{
	this->nlevels = topos.size();
	for (auto i : range<std::size_t>(3)) {
		dims[i] = aarray<int, len_t, 2>(topos[0]->nproc(i), nlevels);
		coord[i] = topos[0]->coord(i);
	}
	line_data = std::make_unique<line_pkg>(*topos[0]);

	init_gfunc(topos);
	init_so(topos);

	for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
		for (int dir = 0; dir < halo_dir::count; dir++) {
			send_active[index(lvl,dir)] = true;
			recv_active[index(lvl,dir)] = true;
		}
		if (not params->periodic[0]) {
			if (topos[0]->coord(0) == 0) {
				send_active[index(lvl, halo_dir::west)] = false;
				recv_active[index(lvl, halo_dir::east)] = false;
			}
			if (topos[0]->coord(0) == topos[0]->nproc(0) - 1) {
				send_active[index(lvl, halo_dir::east)] = false;
				recv_active[index(lvl, halo_dir::west)] = false;
			}
		}
		if (not params->periodic[1]) {
			if (topos[0]->coord(1) == 0) {
				send_active[index(lvl, halo_dir::south)] = false;
				recv_active[index(lvl, halo_dir::north)] = false;
			}
			if (topos[0]->coord(1) == topos[0]->nproc(1) - 1) {
				send_active[index(lvl, halo_dir::north)] = false;
				recv_active[index(lvl, halo_dir::south)] = false;
			}
		}
		if (not params->periodic[2]) {
			if (topos[0]->coord(2) == 0) {
				send_active[index(lvl, halo_dir::bottom)] = false;
				recv_active[index(lvl, halo_dir::top)] = false;
			}
			if (topos[0]->coord(2) == topos[0]->nproc(2) - 1) {
				send_active[index(lvl, halo_dir::top)] = false;
				recv_active[index(lvl, halo_dir::bottom)] = false;
			}
		}
	}

	init_dims(*topos[0]);
}

void tausch_exchanger::init_dims(grid_topo & topo)
{
	auto & dimxfine = dimfine[0];
	auto & dimyfine = dimfine[1];
	auto & dimzfine = dimfine[2];

	if (topo.dimzfine.size() > 0) {
		dimzfine = topo.dimzfine;
	} else {
		dimzfine.reserve(topo.nproc(2));
		for (auto k : range(topo.nproc(2))) {
			dimzfine[k] = topo.nlocal(2) - nghost[2];
		}
	}

	if (topo.dimyfine.size() > 0) {
		dimyfine = topo.dimyfine;
	} else {
		dimyfine.reserve(topo.nproc(1));
		for (auto j : range(topo.nproc(1))) {
			dimyfine[j] = topo.nlocal(1) - nghost[1];
		}
	}

	if (topo.dimxfine.size() > 0) {
		dimxfine = topo.dimxfine;
	} else {
		dimxfine.reserve(topo.nproc(0));
		for (auto i : range(topo.nproc(0))) {
			dimxfine[i] = topo.nlocal(0) - nghost[0];
		}
	}

	BMG3_SymStd_SETUP_Tausch(topo.nlevel(),
	                         dims[0].data(), dims[1].data(), dims[2].data(),
	                         dimxfine.data(), dimyfine.data(), dimzfine.data(),
	                         topo.nproc(0), topo.nproc(1), topo.nproc(2));

	init_datadist();
}


void tausch_exchanger::init_datadist()
{
	auto & datadist = line_data->datadist;

	for (int dim = 0; dim < 3; dim++) {
		datadist[dim](0,0) = 2;
		datadist[dim](1,0) = datadist[dim](0,0) + dimfine[dim][0] - 1;
		for (auto i : range<len_t>(1, datadist[dim].len(1))) {
			datadist[dim](0, i) = datadist[dim](1, i-1) + 1;
			datadist[dim](1, i) = datadist[dim](0, i) + dimfine[dim][i] - 1;
		}
	}

	line_data->linebuf.reserve(std::max({dimfine[0][0], dimfine[1][0], dimfine[2][0]}) * 8 *
	                           std::max({datadist[0].len(1), datadist[1].len(1), datadist[2].len(1)}));
}


void tausch_exchanger::init_gfunc(std::vector<topo_ptr> & topos)
{
	std::vector<TauschHaloSpec> remote_spec;
	std::vector<TauschHaloSpec> local_spec;

	remote_spec.reserve(halo_dir::count * nlevels);
	local_spec.reserve(halo_dir::count * nlevels);
	send_active.reserve(halo_dir::count * nlevels);
	recv_active.reserve(halo_dir::count * nlevels);

	int rank;
	MPI_Comm_rank(topos[0]->comm, &rank);

	for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
		set_level_spec(lvl, rank,
		               *topos[lvl],
		               remote_spec, local_spec);
	}

	int nbuf = 1;

	tausch = std::make_unique<Tausch3D<real_t>>(MPI_DOUBLE, nbuf, nullptr, topos[0]->comm);

	for (std::size_t i = 0; i < halo_dir::count * nlevels; i++) {
		tausch->addLocalHaloInfoCwC(local_spec[i]);
		tausch->addRemoteHaloInfoCwC(remote_spec[i]);
	}
}


void tausch_exchanger::init_so(std::vector<topo_ptr> & topos)
{
	std::vector<TauschHaloSpec> remote_spec;
	std::vector<TauschHaloSpec> local_spec;

	remote_spec.reserve(halo_dir::count * nlevels);
	local_spec.reserve(halo_dir::count * nlevels);
	send_active.reserve(halo_dir::count * nlevels);
	recv_active.reserve(halo_dir::count * nlevels);

	int rank;
	MPI_Comm_rank(topos[0]->comm, &rank);

	for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
		set_level_spec_so(lvl, rank,
		                  *topos[lvl],
		                  remote_spec, local_spec);
	}

	int nbuf = stencil_ndirs<xxvii_pt>::value;

	tausch_so = std::make_unique<Tausch3D<real_t>>(MPI_DOUBLE, nbuf, nullptr, topos[0]->comm);

	for (std::size_t i = 0; i < halo_dir::count * nlevels; i++) {
		tausch_so->addLocalHaloInfoCwC(local_spec[i]);
		tausch_so->addRemoteHaloInfoCwC(remote_spec[i]);
	}
}


void tausch_exchanger::set_level_spec(int lvl, int rank,
                                      grid_topo & topo,
                                      std::vector<TauschHaloSpec> & remote_spec,
                                      std::vector<TauschHaloSpec> & local_spec)
{
	for (int i = 0; i < halo_dir::count; i++) {
		remote_spec[index(lvl,i)].bufferWidth = topo.nlocal(0);
		remote_spec[index(lvl,i)].bufferHeight = topo.nlocal(1);
		remote_spec[index(lvl,i)].bufferDepth = topo.nlocal(2);
		local_spec[index(lvl,i)].bufferWidth = topo.nlocal(0);
		local_spec[index(lvl,i)].bufferHeight = topo.nlocal(1);
		local_spec[index(lvl,i)].bufferDepth = topo.nlocal(2);
	}

	// right
	remote_spec[index(lvl,halo_dir::east)].haloX = 0;
	remote_spec[index(lvl,halo_dir::east)].haloY = 0;
	remote_spec[index(lvl,halo_dir::east)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::east)].haloWidth = 1;
	remote_spec[index(lvl,halo_dir::east)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::east)].haloDepth  = topo.nlocal(2);
	if (topo.coord(0) == 0)
		remote_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank + topo.nproc(0) - 1;
	else
		remote_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank - 1;

	local_spec[index(lvl,halo_dir::east)].haloX = topo.nlocal(0)-2;
	local_spec[index(lvl,halo_dir::east)].haloY = 0;
	local_spec[index(lvl,halo_dir::east)].haloZ = 0;
	local_spec[index(lvl,halo_dir::east)].haloWidth = 1;
	local_spec[index(lvl,halo_dir::east)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::east)].haloDepth  = topo.nlocal(2);
	if (topo.coord(0) == (topo.nproc(0) - 1))
		local_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank - topo.nproc(0) + 1;
	else
		local_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank + 1;


	// left
	remote_spec[index(lvl,halo_dir::west)].haloX = topo.nlocal(0) - 1;
	remote_spec[index(lvl,halo_dir::west)].haloY = 0;
	remote_spec[index(lvl,halo_dir::west)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::west)].haloWidth = 1;
	remote_spec[index(lvl,halo_dir::west)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::west)].haloDepth  = topo.nlocal(2);
	remote_spec[index(lvl,halo_dir::west)].remoteMpiRank = local_spec[index(lvl,halo_dir::east)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::west)].haloX = 1;
	local_spec[index(lvl,halo_dir::west)].haloY = 0;
	local_spec[index(lvl,halo_dir::west)].haloZ = 0;
	local_spec[index(lvl,halo_dir::west)].haloWidth = 1;
	local_spec[index(lvl,halo_dir::west)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::west)].haloDepth  = topo.nlocal(2);
	local_spec[index(lvl,halo_dir::west)].remoteMpiRank = remote_spec[index(lvl,halo_dir::east)].remoteMpiRank;


	// north
	remote_spec[index(lvl,halo_dir::north)].haloX = 0;
	remote_spec[index(lvl,halo_dir::north)].haloY = 0;
	remote_spec[index(lvl,halo_dir::north)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::north)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::north)].haloHeight = 1;
	remote_spec[index(lvl,halo_dir::north)].haloDepth  = topo.nlocal(2);
	if (topo.coord(1) == 0)
		remote_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
	else
		remote_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank - topo.nproc(0);

	local_spec[index(lvl,halo_dir::north)].haloX = 0;
	local_spec[index(lvl,halo_dir::north)].haloY = topo.nlocal(1) - 2;
	local_spec[index(lvl,halo_dir::north)].haloZ = 0;
	local_spec[index(lvl,halo_dir::north)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::north)].haloHeight = 1;
	local_spec[index(lvl,halo_dir::north)].haloDepth  = topo.nlocal(2);
	if (topo.coord(1) == (topo.nproc(1) - 1))
		local_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
	else
		local_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank + topo.nproc(0);
	// south
	remote_spec[index(lvl,halo_dir::south)].haloX = 0;
	remote_spec[index(lvl,halo_dir::south)].haloY = topo.nlocal(1) - 1;
	remote_spec[index(lvl,halo_dir::south)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::south)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::south)].haloHeight = 1;
	remote_spec[index(lvl,halo_dir::south)].haloDepth  = topo.nlocal(2);
	remote_spec[index(lvl,halo_dir::south)].remoteMpiRank = local_spec[index(lvl,halo_dir::north)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::south)].haloX = 0;
	local_spec[index(lvl,halo_dir::south)].haloY = 1;
	local_spec[index(lvl,halo_dir::south)].haloZ = 0;
	local_spec[index(lvl,halo_dir::south)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::south)].haloHeight = 1;
	local_spec[index(lvl,halo_dir::south)].haloDepth  = topo.nlocal(2);
	local_spec[index(lvl,halo_dir::south)].remoteMpiRank = remote_spec[index(lvl,halo_dir::north)].remoteMpiRank;


	// top
	remote_spec[index(lvl,halo_dir::top)].haloX = 0;
	remote_spec[index(lvl,halo_dir::top)].haloY = 0;
	remote_spec[index(lvl,halo_dir::top)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::top)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::top)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::top)].haloDepth  = 1;
	if (topo.coord(2) == 0) {
		remote_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank
			+ topo.nproc(0)*topo.nproc(1)*topo.nproc(2)
			- topo.nproc(0)*topo.nproc(1);
	} else {
		remote_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1);
	}

	local_spec[index(lvl,halo_dir::top)].haloX = 0;
	local_spec[index(lvl,halo_dir::top)].haloY = 0;
	local_spec[index(lvl,halo_dir::top)].haloZ = topo.nlocal(2) - 2;
	local_spec[index(lvl,halo_dir::top)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::top)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::top)].haloDepth  = 1;
	if (topo.coord(2) == (topo.nproc(2) - 1)) {
		local_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank
			- topo.nproc(0)*topo.nproc(1)*topo.nproc(2)
			+ topo.nproc(0)*topo.nproc(1);
	} else {
		local_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1);
	}


	// bottom
	remote_spec[index(lvl,halo_dir::bottom)].haloX = 0;
	remote_spec[index(lvl,halo_dir::bottom)].haloY = 0;
	remote_spec[index(lvl,halo_dir::bottom)].haloZ = topo.nlocal(2) - 1;
	remote_spec[index(lvl,halo_dir::bottom)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::bottom)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::bottom)].haloDepth  = 1;
	remote_spec[index(lvl,halo_dir::bottom)].remoteMpiRank = local_spec[index(lvl,halo_dir::top)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::bottom)].haloX = 0;
	local_spec[index(lvl,halo_dir::bottom)].haloY = 0;
	local_spec[index(lvl,halo_dir::bottom)].haloZ = 1;
	local_spec[index(lvl,halo_dir::bottom)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::bottom)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::bottom)].haloDepth  = 1;
	local_spec[index(lvl,halo_dir::bottom)].remoteMpiRank = remote_spec[index(lvl,halo_dir::top)].remoteMpiRank;
}


void tausch_exchanger::set_level_spec_so(int lvl, int rank,
                                         grid_topo & topo,
                                         std::vector<TauschHaloSpec> & remote_spec,
                                         std::vector<TauschHaloSpec> & local_spec)
{
	for (int i = 0; i < halo_dir::count; i++) {
		remote_spec[index(lvl,i)].bufferWidth = topo.nlocal(0) + 1;
		remote_spec[index(lvl,i)].bufferHeight = topo.nlocal(1) + 1;
		remote_spec[index(lvl,i)].bufferDepth = topo.nlocal(2) + 1;
		local_spec[index(lvl,i)].bufferWidth = topo.nlocal(0) + 1;
		local_spec[index(lvl,i)].bufferHeight = topo.nlocal(1) + 1;
		local_spec[index(lvl,i)].bufferDepth = topo.nlocal(2) + 1;
	}

	// right
	remote_spec[index(lvl,halo_dir::east)].haloX = 0;
	remote_spec[index(lvl,halo_dir::east)].haloY = 0;
	remote_spec[index(lvl,halo_dir::east)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::east)].haloWidth = 1;
	remote_spec[index(lvl,halo_dir::east)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::east)].haloDepth  = topo.nlocal(2);
	if (topo.coord(0) == 0)
		remote_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank + topo.nproc(0) - 1;
	else
		remote_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank - 1;

	local_spec[index(lvl,halo_dir::east)].haloX = topo.nlocal(0)-2;
	local_spec[index(lvl,halo_dir::east)].haloY = 0;
	local_spec[index(lvl,halo_dir::east)].haloZ = 0;
	local_spec[index(lvl,halo_dir::east)].haloWidth = 1;
	local_spec[index(lvl,halo_dir::east)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::east)].haloDepth  = topo.nlocal(2);
	if (topo.coord(0) == (topo.nproc(0) - 1))
		local_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank - topo.nproc(0) + 1;
	else
		local_spec[index(lvl,halo_dir::east)].remoteMpiRank = rank + 1;


	// left
	remote_spec[index(lvl,halo_dir::west)].haloX = topo.nlocal(0) - 1;
	remote_spec[index(lvl,halo_dir::west)].haloY = 0;
	remote_spec[index(lvl,halo_dir::west)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::west)].haloWidth = 2;
	remote_spec[index(lvl,halo_dir::west)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::west)].haloDepth  = topo.nlocal(2);
	remote_spec[index(lvl,halo_dir::west)].remoteMpiRank = local_spec[index(lvl,halo_dir::east)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::west)].haloX = 1;
	local_spec[index(lvl,halo_dir::west)].haloY = 0;
	local_spec[index(lvl,halo_dir::west)].haloZ = 0;
	local_spec[index(lvl,halo_dir::west)].haloWidth = 2;
	local_spec[index(lvl,halo_dir::west)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::west)].haloDepth  = topo.nlocal(2);
	local_spec[index(lvl,halo_dir::west)].remoteMpiRank = remote_spec[index(lvl,halo_dir::east)].remoteMpiRank;


	// north
	remote_spec[index(lvl,halo_dir::north)].haloX = 0;
	remote_spec[index(lvl,halo_dir::north)].haloY = 0;
	remote_spec[index(lvl,halo_dir::north)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::north)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::north)].haloHeight = 1;
	remote_spec[index(lvl,halo_dir::north)].haloDepth  = topo.nlocal(2);
	if (topo.coord(1) == 0)
		remote_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
	else
		remote_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank - topo.nproc(0);

	local_spec[index(lvl,halo_dir::north)].haloX = 0;
	local_spec[index(lvl,halo_dir::north)].haloY = topo.nlocal(1) - 2;
	local_spec[index(lvl,halo_dir::north)].haloZ = 0;
	local_spec[index(lvl,halo_dir::north)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::north)].haloHeight = 1;
	local_spec[index(lvl,halo_dir::north)].haloDepth  = topo.nlocal(2);
	if (topo.coord(1) == (topo.nproc(1) - 1))
		local_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
	else
		local_spec[index(lvl,halo_dir::north)].remoteMpiRank = rank + topo.nproc(0);
	// south
	remote_spec[index(lvl,halo_dir::south)].haloX = 0;
	remote_spec[index(lvl,halo_dir::south)].haloY = topo.nlocal(1) - 1;
	remote_spec[index(lvl,halo_dir::south)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::south)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::south)].haloHeight = 2;
	remote_spec[index(lvl,halo_dir::south)].haloDepth  = topo.nlocal(2);
	remote_spec[index(lvl,halo_dir::south)].remoteMpiRank = local_spec[index(lvl,halo_dir::north)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::south)].haloX = 0;
	local_spec[index(lvl,halo_dir::south)].haloY = 1;
	local_spec[index(lvl,halo_dir::south)].haloZ = 0;
	local_spec[index(lvl,halo_dir::south)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::south)].haloHeight = 2;
	local_spec[index(lvl,halo_dir::south)].haloDepth  = topo.nlocal(2);
	local_spec[index(lvl,halo_dir::south)].remoteMpiRank = remote_spec[index(lvl,halo_dir::north)].remoteMpiRank;


	// top
	remote_spec[index(lvl,halo_dir::top)].haloX = 0;
	remote_spec[index(lvl,halo_dir::top)].haloY = 0;
	remote_spec[index(lvl,halo_dir::top)].haloZ = 0;
	remote_spec[index(lvl,halo_dir::top)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::top)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::top)].haloDepth  = 1;
	if (topo.coord(2) == 0) {
		remote_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank
			+ topo.nproc(0)*topo.nproc(1)*topo.nproc(2)
			- topo.nproc(0)*topo.nproc(1);
	} else {
		remote_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1);
	}

	local_spec[index(lvl,halo_dir::top)].haloX = 0;
	local_spec[index(lvl,halo_dir::top)].haloY = 0;
	local_spec[index(lvl,halo_dir::top)].haloZ = topo.nlocal(2) - 2;
	local_spec[index(lvl,halo_dir::top)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::top)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::top)].haloDepth  = 1;
	if (topo.coord(2) == (topo.nproc(2) - 1)) {
		local_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank
			- topo.nproc(0)*topo.nproc(1)*topo.nproc(2)
			+ topo.nproc(0)*topo.nproc(1);
	} else {
		local_spec[index(lvl,halo_dir::top)].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1);
	}


	// bottom
	remote_spec[index(lvl,halo_dir::bottom)].haloX = 0;
	remote_spec[index(lvl,halo_dir::bottom)].haloY = 0;
	remote_spec[index(lvl,halo_dir::bottom)].haloZ = topo.nlocal(2) - 1;
	remote_spec[index(lvl,halo_dir::bottom)].haloWidth  = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::bottom)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::bottom)].haloDepth  = 2;
	remote_spec[index(lvl,halo_dir::bottom)].remoteMpiRank = local_spec[index(lvl,halo_dir::top)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::bottom)].haloX = 0;
	local_spec[index(lvl,halo_dir::bottom)].haloY = 0;
	local_spec[index(lvl,halo_dir::bottom)].haloZ = 1;
	local_spec[index(lvl,halo_dir::bottom)].haloWidth  = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::bottom)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::bottom)].haloDepth  = 2;
	local_spec[index(lvl,halo_dir::bottom)].remoteMpiRank = remote_spec[index(lvl,halo_dir::top)].remoteMpiRank;
}


void tausch_exchanger::exchange_func(int k, real_t * gf)
{
	int lvl = nlevels - k;
	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (recv_active[index(lvl, dir)])
			tausch->postReceiveCwC(index(lvl, dir), index(lvl, dir));
	}

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)]) {
			tausch->packSendBufferCwC(index(lvl,dir), 0, gf);
			tausch->sendCwC(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch->recvCwC(index(lvl,dir));
			tausch->unpackRecvBufferCwC(index(lvl,dir), 0, gf);
		}
	}
}

void tausch_exchanger::exchange_sten(int k, real_t * so)
{
	int lvl = nlevels - k;

	auto & dimx = leveldims(0);
	auto & dimy = leveldims(1);
	auto & dimz = leveldims(2);

	len_t II = dimx(coord[0], k-1) + 2;
	len_t JJ = dimy(coord[1], k-1) + 2;
	len_t KK = dimz(coord[2], k-1) + 2;

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (recv_active[index(lvl, dir)])
			tausch_so->postReceiveCwC(index(lvl, dir), index(lvl, dir));
	}

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)]) {
			for (int sdir = 0; sdir < stencil_ndirs<xxvii_pt>::value; sdir++) {
				len_t offset = (JJ+1)*(II+1)*(KK+1)*sdir;
				tausch_so->packSendBufferCwC(index(lvl,dir), sdir, so + offset);
			}
			tausch_so->sendCwC(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch_so->recvCwC(index(lvl,dir));
			for (int sdir = 0; sdir < stencil_ndirs<xxvii_pt>::value; sdir++) {
				len_t offset = (JJ+1)*(II+1)*(KK+1)*sdir;
				tausch_so->unpackRecvBufferCwC(index(lvl,dir), sdir, so + offset);
			}
		}
	}
}

void tausch_exchanger::run(mpi::grid_func & f)
{
	auto lvl = f.grid().level();
	lvl = nlevels - lvl - 1;

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (recv_active[index(lvl, dir)])
			tausch->postReceiveCwC(index(lvl, dir), index(lvl, dir));
	}

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)]) {
			tausch->packSendBufferCwC(index(lvl,dir), 0, f.data());
			tausch->sendCwC(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch->recvCwC(index(lvl,dir));
			tausch->unpackRecvBufferCwC(index(lvl,dir), 0, f.data());
		}
	}
}

void tausch_exchanger::activate_send(halo_dir dir, bool active)
{
	for (auto lvl : range<std::size_t>(nlevels))
		send_active[index(lvl, dir)] = active;
}


void tausch_exchanger::activate_recv(halo_dir dir, bool active)
{
	for (auto lvl : range<std::size_t>(nlevels))
		recv_active[index(lvl, dir)] = active;
}

}}}
