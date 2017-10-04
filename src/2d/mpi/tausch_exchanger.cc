#include <cedar/2d/mpi/tausch_exchanger.h>


namespace cedar { namespace cdr2 { namespace mpi {

tausch_exchanger::tausch_exchanger(const kernel_params & params,
                                   std::vector<topo_ptr> topos) : nlevels(topos.size())
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
		for (int dir = 0; dir < halo_dir::count; dir++)  {
			send_active[index(lvl,dir)] = true;
			recv_active[index(lvl,dir)] = true;
		}
		if (not params.periodic[0]) {
			if (topos[0]->coord(0) == 0) {
				send_active[index(lvl, halo_dir::left)] = false;
				recv_active[index(lvl, halo_dir::right)] = false;
			}
			if (topos[0]->coord(0) == topos[0]->nproc(0) - 1) {
				send_active[index(lvl, halo_dir::right)] = false;
				recv_active[index(lvl, halo_dir::left)] = false;
			}
		}
		if (not params.periodic[1]) {
			if (topos[0]->coord(1) == 0) {
				send_active[index(lvl, halo_dir::down)] = false;
				recv_active[index(lvl, halo_dir::up)] = false;
			}
			if (topos[0]->coord(1) == topos[0]->nproc(1) - 1) {
				send_active[index(lvl, halo_dir::up)] = false;
				recv_active[index(lvl, halo_dir::down)] = false;
			}
		}
	}

	int nbuf = 1;

	tausch = std::make_unique<Tausch2D<real_t>>(MPI_DOUBLE, nbuf, nullptr, topos[0]->comm);

	tausch->setLocalHaloInfo(TAUSCH_CwC, halo_dir::count * nlevels, local_spec.data());
	tausch->setRemoteHaloInfo(TAUSCH_CwC, halo_dir::count * nlevels, remote_spec.data());
}


void tausch_exchanger::set_level_spec(int lvl, int rank,
                                      grid_topo & topo,
                                      std::vector<TauschHaloSpec> & remote_spec,
                                      std::vector<TauschHaloSpec> & local_spec)
{
	for (int i = 0; i < halo_dir::count; i++) {
		remote_spec[index(lvl,i)].bufferWidth = topo.nlocal(0);
		remote_spec[index(lvl,i)].bufferHeight = topo.nlocal(1);
		local_spec[index(lvl,i)].bufferWidth = topo.nlocal(0);
		local_spec[index(lvl,i)].bufferHeight = topo.nlocal(1);
	}

	// right
	remote_spec[index(lvl,halo_dir::right)].haloX = 0;
	remote_spec[index(lvl,halo_dir::right)].haloY = 0;
	remote_spec[index(lvl,halo_dir::right)].haloWidth = 1;
	remote_spec[index(lvl,halo_dir::right)].haloHeight = topo.nlocal(1);
	if (topo.coord(0) == 0)
		remote_spec[index(lvl,halo_dir::right)].remoteMpiRank = rank + topo.nproc(0) - 1;
	else
		remote_spec[index(lvl,halo_dir::right)].remoteMpiRank = rank - 1;

	local_spec[index(lvl,halo_dir::right)].haloX = topo.nlocal(0)-2;
	local_spec[index(lvl,halo_dir::right)].haloY = 0;
	local_spec[index(lvl,halo_dir::right)].haloWidth = 1;
	local_spec[index(lvl,halo_dir::right)].haloHeight = topo.nlocal(1);
	if ((rank + 1) % topo.nproc(0) == 0)
		local_spec[index(lvl,halo_dir::right)].remoteMpiRank = rank - topo.nproc(0) + 1;
	else
		local_spec[index(lvl,halo_dir::right)].remoteMpiRank = rank + 1;


	// left
	remote_spec[index(lvl,halo_dir::left)].haloX = topo.nlocal(0) - 1;
	remote_spec[index(lvl,halo_dir::left)].haloY = 0;
	remote_spec[index(lvl,halo_dir::left)].haloWidth = 1;
	remote_spec[index(lvl,halo_dir::left)].haloHeight = topo.nlocal(1);
	remote_spec[index(lvl,halo_dir::left)].remoteMpiRank = local_spec[index(lvl,halo_dir::right)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::left)].haloX = 1;
	local_spec[index(lvl,halo_dir::left)].haloY = 0;
	local_spec[index(lvl,halo_dir::left)].haloWidth = 1;
	local_spec[index(lvl,halo_dir::left)].haloHeight = topo.nlocal(1);
	local_spec[index(lvl,halo_dir::left)].remoteMpiRank = remote_spec[index(lvl,halo_dir::right)].remoteMpiRank;


	// up
	remote_spec[index(lvl,halo_dir::up)].haloX = 0;
	remote_spec[index(lvl,halo_dir::up)].haloY = 0;
	remote_spec[index(lvl,halo_dir::up)].haloWidth = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::up)].haloHeight = 1;
	if (topo.coord(1) == 0)
		remote_spec[index(lvl,halo_dir::up)].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
	else
		remote_spec[index(lvl,halo_dir::up)].remoteMpiRank = rank - topo.nproc(0);

	local_spec[index(lvl,halo_dir::up)].haloX = 0;
	local_spec[index(lvl,halo_dir::up)].haloY = topo.nproc(1) - 2;
	local_spec[index(lvl,halo_dir::up)].haloWidth = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::up)].haloHeight = 1;
	if (topo.coord(1) == (topo.nproc(1) - 1))
		local_spec[index(lvl,halo_dir::up)].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
	else
		local_spec[index(lvl,halo_dir::up)].remoteMpiRank = rank + topo.nproc(0);
	// down
	remote_spec[index(lvl,halo_dir::down)].haloX = 0;
	remote_spec[index(lvl,halo_dir::down)].haloY = topo.nlocal(1) - 1;
	remote_spec[index(lvl,halo_dir::down)].haloWidth = topo.nlocal(0);
	remote_spec[index(lvl,halo_dir::down)].haloHeight = 1;
	remote_spec[index(lvl,halo_dir::down)].remoteMpiRank = local_spec[index(lvl,halo_dir::up)].remoteMpiRank;

	local_spec[index(lvl,halo_dir::down)].haloX = 0;
	local_spec[index(lvl,halo_dir::down)].haloY = 1;
	local_spec[index(lvl,halo_dir::down)].haloWidth = topo.nlocal(0);
	local_spec[index(lvl,halo_dir::down)].haloHeight = 1;
	local_spec[index(lvl,halo_dir::down)].remoteMpiRank = remote_spec[index(lvl,halo_dir::up)].remoteMpiRank;
}


void tausch_exchanger::exchange_func(int k, real_t * gf)
{
}

void tausch_exchanger::exchange_sten(int k, real_t * so)
{
}

void tausch_exchanger::exchange(mpi::grid_func & f)
{
	auto lvl = f.grid().level();

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (recv_active[index(lvl, dir)])
			tausch->postReceive(TAUSCH_CwC, index(lvl, dir), index(lvl, dir));
	}

	for (int dir = 0; dir < halo_dir::count; dir++) {
		if (send_active[index(lvl, dir)]) {
			tausch->packSendBuffer(TAUSCH_CwC, index(lvl,dir), 0, f.data());
			tausch->send(TAUSCH_CwC, index(lvl,dir), index(lvl,dir));

		}
		if (recv_active[index(lvl, dir)]) {
			tausch->recv(TAUSCH_CwC, index(lvl,dir));
			tausch->unpackRecvBuffer(TAUSCH_CwC, index(lvl,dir), 0, f.data());
		}
	}
}


}}}
