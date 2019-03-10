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


void tausch_exchanger::setup_worker(std::vector<topo_ptr> topos)
{
	this->nlevels = topos.size();
	for (auto i : range<std::size_t>(3)) {
		dims[i] = aarray<int, len_t, 2>(topos[0]->nproc(i), nlevels);
		coord[i] = topos[0]->coord(i);
	}
	line_data = std::make_unique<line_pkg>(*topos[0]);

	init_dims(*topos[0]);
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
	std::vector<std::vector<std::array<int, 3> > > remote_spec(halo_dir::count * nlevels);
	std::vector<std::vector<std::array<int, 3> > > local_spec(halo_dir::count * nlevels);
    std::vector<int> remote_remoteMpiRank(halo_dir::count * nlevels);
    std::vector<int> local_remoteMpiRank(halo_dir::count * nlevels);

	send_active.reserve(halo_dir::count * nlevels);
	recv_active.reserve(halo_dir::count * nlevels);

	int rank;
	MPI_Comm_rank(topos[0]->comm, &rank);

	for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
		set_level_spec(lvl, rank,
		               *topos[lvl],
		               remote_spec, local_spec,
                       remote_remoteMpiRank, local_remoteMpiRank);
	}

	int nbuf = 1;

	tausch = std::make_unique<Tausch<real_t>>(MPI_DOUBLE, topos[0]->comm);

	for (std::size_t i = 0; i < halo_dir::count * nlevels; i++) {
		tausch->addLocalHaloInfo(local_spec[i], nbuf, local_remoteMpiRank[i]);
		tausch->addRemoteHaloInfo(remote_spec[i], nbuf, remote_remoteMpiRank[i]);
	}
}


void tausch_exchanger::init_so(std::vector<topo_ptr> & topos)
{
	std::vector<std::vector<std::array<int, 3> > > remote_spec(halo_dir::count * nlevels);
	std::vector<std::vector<std::array<int, 3> > > local_spec(halo_dir::count * nlevels);
    std::vector<int> remote_remoteMpiRank(halo_dir::count * nlevels);
    std::vector<int> local_remoteMpiRank(halo_dir::count * nlevels);

	send_active.reserve(halo_dir::count * nlevels);
	recv_active.reserve(halo_dir::count * nlevels);

	int rank;
	MPI_Comm_rank(topos[0]->comm, &rank);

	for (std::size_t lvl = 0; lvl < nlevels; lvl++) {
		set_level_spec_so(lvl, rank,
		                  *topos[lvl],
		                  remote_spec, local_spec,
                          remote_remoteMpiRank, local_remoteMpiRank);
	}

	int nbuf = stencil_ndirs<xxvii_pt>::value;

	tausch_so = std::make_unique<Tausch<real_t>>(MPI_DOUBLE, topos[0]->comm);

	for (std::size_t i = 0; i < halo_dir::count * nlevels; i++) {
		tausch_so->addLocalHaloInfo(local_spec[i], nbuf, local_remoteMpiRank[i]);
		tausch_so->addRemoteHaloInfo(remote_spec[i], nbuf, remote_remoteMpiRank[i]);
	}
}


void tausch_exchanger::set_level_spec(int lvl, int rank,
                                      grid_topo & topo,
                                      std::vector<std::vector<std::array<int, 3> > > & remote_spec,
                                      std::vector<std::vector<std::array<int, 3> > > & local_spec,
                                      std::vector<int> & remote_remoteMpiRank,
                                      std::vector<int> & local_remoteMpiRank)
{
	const int nx = topo.nlocal(0);
    const int ny = topo.nlocal(1);
    const int nz = topo.nlocal(2);

	// right

    for(int z = 0; z < nz; ++z) {
        remote_spec[index(lvl,halo_dir::east)].push_back(std::array<int, 3>{z*(nx*ny)       , ny, nx});
        local_spec [index(lvl,halo_dir::east)].push_back(std::array<int, 3>{z*(nx*ny) + nx-2, ny, nx});
    }

    if(topo.coord(0) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::east)] = rank + topo.nproc(0)-1;
    else
        remote_remoteMpiRank[index(lvl,halo_dir::east)] = rank - 1;

    if(topo.coord(0) == (topo.nproc(0) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::east)] = rank - topo.nproc(0) + 1;
    else
        local_remoteMpiRank[index(lvl,halo_dir::east)] = rank + 1;

	// left

    for(int z = 0; z < nz; ++z) {
        remote_spec[index(lvl,halo_dir::west)].push_back(std::array<int, 3>{z*(nx*ny) + nx-1, ny, nx});
        local_spec [index(lvl,halo_dir::west)].push_back(std::array<int, 3>{z*(nx*ny) + 1   , ny, nx});
    }

    remote_remoteMpiRank[index(lvl,halo_dir::west)] = local_remoteMpiRank[index(lvl,halo_dir::east)];
    local_remoteMpiRank[index(lvl,halo_dir::west)] = remote_remoteMpiRank[index(lvl,halo_dir::east)];

	// north

    for(int z = 0; z < nz; ++z) {
        remote_spec[index(lvl,halo_dir::north)].push_back(std::array<int, 3>{z*(nx*ny)            , nx, 1});
        local_spec [index(lvl,halo_dir::north)].push_back(std::array<int, 3>{z*(nx*ny) + (ny-2)*nx, nx, 1});
    }

    if(topo.coord(1) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::north)] = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
    else
        remote_remoteMpiRank[index(lvl,halo_dir::north)] = rank - topo.nproc(0);

    if(topo.coord(1) == (topo.nproc(1) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::north)] = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
    else
        local_remoteMpiRank[index(lvl,halo_dir::north)] = rank + topo.nproc(0);

	// south

    for(int z = 0; z < nz; ++z) {
        remote_spec[index(lvl,halo_dir::south)].push_back(std::array<int, 3>{z*(nx*ny) + (ny-1)*nx, nx, 1});
        local_spec [index(lvl,halo_dir::south)].push_back(std::array<int, 3>{z*(nx*ny) +        nx, nx, 1});
    }

    remote_remoteMpiRank[index(lvl,halo_dir::south)] = local_remoteMpiRank[index(lvl,halo_dir::north)];
    local_remoteMpiRank[index(lvl,halo_dir::south)] = remote_remoteMpiRank[index(lvl,halo_dir::north)];

	// top

    for(int y = 0; y < ny; ++y) {
        remote_spec[index(lvl,halo_dir::top)].push_back(std::array<int, 3>{                 y*nx, nx, 1});
        local_spec [index(lvl,halo_dir::top)].push_back(std::array<int, 3>{(nz-2)*(nx*ny) + y*nx, nx, 1});
    }

    if(topo.coord(2) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::top)] = rank + topo.nproc(0)*topo.nproc(1)*topo.nproc(2) - topo.nproc(0)*topo.nproc(1);
    else
        remote_remoteMpiRank[index(lvl,halo_dir::top)] = rank - topo.nproc(0)*topo.nproc(1);

    if(topo.coord(2) == (topo.nproc(2) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::top)] = rank - topo.nproc(0)*topo.nproc(1)*topo.nproc(2) + topo.nproc(0)*topo.nproc(1);
    else
        local_remoteMpiRank[index(lvl,halo_dir::top)] = rank + topo.nproc(0)*topo.nproc(1);

	// bottom

    for(int y = 0; y < ny; ++y) {
        remote_spec[index(lvl,halo_dir::bottom)].push_back(std::array<int, 3>{(nz-1)*(nx*ny) + y*nx, nx, 1});
        local_spec [index(lvl,halo_dir::bottom)].push_back(std::array<int, 3>{       (nx*ny) + y*nx, nx, 1});
    }

    remote_remoteMpiRank[index(lvl,halo_dir::bottom)] = local_remoteMpiRank[index(lvl,halo_dir::top)];
    local_remoteMpiRank[index(lvl,halo_dir::bottom)] = remote_remoteMpiRank[index(lvl,halo_dir::top)];
}


void tausch_exchanger::set_level_spec_so(int lvl, int rank,
                                         grid_topo & topo,
                                         std::vector<std::vector<std::array<int,3> > > & remote_spec,
                                         std::vector<std::vector<std::array<int,3> > > & local_spec,
                                         std::vector<int> & remote_remoteMpiRank,
                                         std::vector<int> & local_remoteMpiRank)
{
	const int nx = topo.nlocal(0)+1;
    const int ny = topo.nlocal(1)+1;
    const int nz = topo.nlocal(2)+1;

	// right

    for(int z = 0; z < nz-1; ++z) {
        remote_spec[index(lvl,halo_dir::east)].push_back(std::array<int, 3>{z*(nx*ny)       , ny-1, nx});
        local_spec [index(lvl,halo_dir::east)].push_back(std::array<int, 3>{z*(nx*ny) + nx-3, ny-1, nx});
    }

    if(topo.coord(0) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::east)] = rank + topo.nproc(0)-1;
    else
        remote_remoteMpiRank[index(lvl,halo_dir::east)] = rank - 1;

    if(topo.coord(0) == (topo.nproc(0) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::east)] = rank - topo.nproc(0) + 1;
    else
        local_remoteMpiRank[index(lvl,halo_dir::east)] = rank + 1;

	// left

    for(int z = 0; z < nz-1; ++z) {
        for(int y = 0; y < ny-1; ++y) {
            remote_spec[index(lvl,halo_dir::west)].push_back(std::array<int, 3>{z*(nx*ny) + y*nx + nx-2, 2, 1});
            local_spec [index(lvl,halo_dir::west)].push_back(std::array<int, 3>{z*(nx*ny) + y*nx + 1   , 2, 1});
        }
    }

    remote_remoteMpiRank[index(lvl,halo_dir::west)] = local_remoteMpiRank[index(lvl,halo_dir::east)];

    // north

    for(int z = 0; z < nz-1; ++z) {
        remote_spec[index(lvl,halo_dir::north)].push_back(std::array<int, 3>{z*(nx*ny)            , nx-1, 1});
        local_spec [index(lvl,halo_dir::north)].push_back(std::array<int, 3>{z*(nx*ny) + (ny-3)*nx, nx-1, 1});
    }

    if(topo.coord(1) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::north)] = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
    else
        remote_remoteMpiRank[index(lvl,halo_dir::north)] = rank - topo.nproc(0);

    if(topo.coord(1) == (topo.nproc(1) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::north)] = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
    else
        local_remoteMpiRank[index(lvl,halo_dir::north)] = rank + topo.nproc(0);

    // south

    for(int z = 0; z < nz-1; ++z) {
        for(int y = 0; y < 2; ++y) {
            remote_spec[index(lvl,halo_dir::south)].push_back(std::array<int, 3>{z*(nx*ny) + (ny-2+y)*nx, nx-1, 1});
            local_spec [index(lvl,halo_dir::south)].push_back(std::array<int, 3>{z*(nx*ny) +    (1+y)*nx, nx-1, 1});
        }
    }

    remote_remoteMpiRank[index(lvl,halo_dir::south)] = local_remoteMpiRank[index(lvl,halo_dir::north)];
    local_remoteMpiRank[index(lvl,halo_dir::south)] = remote_remoteMpiRank[index(lvl,halo_dir::north)];

    // top

    for(int y = 0; y < ny-1; ++y) {
        remote_spec[index(lvl,halo_dir::top)].push_back(std::array<int, 3>{                 y*nx, nx-1, 1});
        local_spec [index(lvl,halo_dir::top)].push_back(std::array<int, 3>{(nz-3)*(nx*ny) + y*nx, nx-1, 1});
    }

    if(topo.coord(2) == 0)
        remote_remoteMpiRank[index(lvl,halo_dir::top)] = rank + topo.nproc(0)*topo.nproc(1)*topo.nproc(2) - topo.nproc(0)*topo.nproc(1);
    else
        remote_remoteMpiRank[index(lvl,halo_dir::top)] = rank - topo.nproc(0)*topo.nproc(1);

    if(topo.coord(2) == (topo.nproc(2) - 1))
        local_remoteMpiRank[index(lvl,halo_dir::top)] = rank - topo.nproc(0)*topo.nproc(1)*topo.nproc(2) + topo.nproc(0)*topo.nproc(1);
    else
        local_remoteMpiRank[index(lvl,halo_dir::top)] = rank + topo.nproc(0)*topo.nproc(1);


    // bottom

    for(int z = 0; z < 2; ++z) {
        for(int y = 0; y < ny-1; ++y) {
            remote_spec[index(lvl,halo_dir::bottom)].push_back(std::array<int, 3>{(nz-2+z)*(nx*ny) + y*nx, nx-1, 1});
            local_spec [index(lvl,halo_dir::bottom)].push_back(std::array<int, 3>{   (1+z)*(nx*ny) + y*nx, nx-1, 1});
        }
    }

    remote_remoteMpiRank[index(lvl,halo_dir::bottom)] = local_remoteMpiRank[index(lvl,halo_dir::top)];
    local_remoteMpiRank[index(lvl,halo_dir::bottom)] = remote_remoteMpiRank[index(lvl,halo_dir::top)];

}


void tausch_exchanger::exchange_func(int k, real_t * gf)
{
	int lvl = nlevels - k;

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)]) {
			tausch->packSendBuffer(index(lvl,dir), 0, gf);
			tausch->send(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch->recv(index(lvl,dir), index(lvl,dir));
			tausch->unpackRecvBuffer(index(lvl,dir), 0, gf);
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
		if (send_active[index(lvl, dir)]) {
			for (int sdir = 0; sdir < stencil_ndirs<xxvii_pt>::value; sdir++) {
				len_t offset = (JJ+1)*(II+1)*(KK+1)*sdir;
				tausch_so->packSendBuffer(index(lvl,dir), sdir, so + offset);
			}
			tausch_so->send(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch_so->recv(index(lvl,dir), index(lvl,dir));
			for (int sdir = 0; sdir < stencil_ndirs<xxvii_pt>::value; sdir++) {
				len_t offset = (JJ+1)*(II+1)*(KK+1)*sdir;
				tausch_so->unpackRecvBuffer(index(lvl,dir), sdir, so + offset);
			}
		}
	}
}

void tausch_exchanger::run(mpi::grid_func & f, unsigned short dmask)
{
	auto lvl = f.grid().level();
	lvl = nlevels - lvl - 1;

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)] and (dmask & (1 << (dir / 2)))) {
			tausch->packSendBuffer(index(lvl,dir), 0, f.data());
			tausch->send(index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)] and (dmask & (1 << (dir / 2)))) {
			tausch->recv(index(lvl,dir), index(lvl,dir));
			tausch->unpackRecvBuffer(index(lvl,dir), 0, f.data());
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
