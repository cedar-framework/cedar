#include <cedar/2d/mpi/tausch_exchanger.h>


namespace cedar { namespace cdr2 { namespace mpi {

tausch_exchanger::tausch_exchanger(const kernel_params & params,
                                   grid_topo & topo)
{
	std::array<TauschHaloSpec, halo_dir::count> remote_spec;
	std::array<TauschHaloSpec, halo_dir::count> local_spec;

	int rank;
	MPI_Comm_rank(topo.comm, &rank);

	for (int i = 0; i < halo_dir::count; i++) {
		remote_spec[i].bufferWidth = topo.nlocal(0);
		remote_spec[i].bufferHeight = topo.nlocal(1);
		local_spec[i].bufferWidth = topo.nlocal(0);
		local_spec[i].bufferHeight = topo.nlocal(1);
	}

	// right
	remote_spec[halo_dir::right].haloX = 0;
	remote_spec[halo_dir::right].haloY = 0;
	remote_spec[halo_dir::right].haloWidth = 1;
	remote_spec[halo_dir::right].haloHeight = topo.nlocal(1);
	if (topo.coord(0) == 0)
		remote_spec[halo_dir::right].remoteMpiRank = rank + topo.nproc(0) - 1;
	else
		remote_spec[halo_dir::right].remoteMpiRank = rank - 1;

	local_spec[halo_dir::right].haloX = topo.nlocal(0)-2;
	local_spec[halo_dir::right].haloY = 0;
	local_spec[halo_dir::right].haloWidth = 1;
	local_spec[halo_dir::right].haloHeight = topo.nlocal(1);
	if ((rank + 1) % topo.nproc(0) == 0)
		local_spec[halo_dir::right].remoteMpiRank = rank - topo.nproc(0) + 1;
	else
		local_spec[halo_dir::right].remoteMpiRank = rank + 1;


	// left
	remote_spec[halo_dir::left].haloX = topo.nlocal(0) - 1;
	remote_spec[halo_dir::left].haloY = 0;
	remote_spec[halo_dir::left].haloWidth = 1;
	remote_spec[halo_dir::left].haloHeight = topo.nlocal(1);
	remote_spec[halo_dir::left].remoteMpiRank = local_spec[halo_dir::right].remoteMpiRank;

	local_spec[halo_dir::left].haloX = 1;
	local_spec[halo_dir::left].haloY = 0;
	local_spec[halo_dir::left].haloWidth = 1;
	local_spec[halo_dir::left].haloHeight = topo.nlocal(1);
	local_spec[halo_dir::left].remoteMpiRank = remote_spec[halo_dir::right].remoteMpiRank;


	// up
	remote_spec[halo_dir::up].haloX = 0;
	remote_spec[halo_dir::up].haloY = 0;
	remote_spec[halo_dir::up].haloWidth = topo.nlocal(0);
	remote_spec[halo_dir::up].haloHeight = 1;
	if (topo.coord(1) == 0)
		remote_spec[halo_dir::up].remoteMpiRank = rank + topo.nproc(0)*topo.nproc(1) - topo.nproc(0);
	else
		remote_spec[halo_dir::up].remoteMpiRank = rank - topo.nproc(0);

	local_spec[halo_dir::up].haloX = 0;
	local_spec[halo_dir::up].haloY = topo.nproc(1) - 2;
	local_spec[halo_dir::up].haloWidth = topo.nlocal(0);
	local_spec[halo_dir::up].haloHeight = 1;
	if (topo.coord(1) == (topo.nproc(1) - 1))
		local_spec[halo_dir::up].remoteMpiRank = rank - topo.nproc(0)*topo.nproc(1) + topo.nproc(0);
	else
		local_spec[halo_dir::up].remoteMpiRank = rank + topo.nproc(0);
	// down
	remote_spec[halo_dir::down].haloX = 0;
	remote_spec[halo_dir::down].haloY = topo.nlocal(1) - 1;
	remote_spec[halo_dir::down].haloWidth = topo.nlocal(0);
	remote_spec[halo_dir::down].haloHeight = 1;
	remote_spec[halo_dir::down].remoteMpiRank = local_spec[halo_dir::up].remoteMpiRank;

	local_spec[halo_dir::down].haloX = 0;
	local_spec[halo_dir::down].haloY = 1;
	local_spec[halo_dir::down].haloWidth = topo.nlocal(0);
	local_spec[halo_dir::down].haloHeight = 1;
	local_spec[halo_dir::down].remoteMpiRank = remote_spec[halo_dir::up].remoteMpiRank;

	int nbuf = 1;

	tausch = std::make_unique<Tausch2D<real_t>>(MPI_DOUBLE, nbuf, nullptr, topo.comm);

	tausch->setLocalHaloInfo(TAUSCH_CwC, halo_dir::count, local_spec.data());
	tausch->setRemoteHaloInfo(TAUSCH_CwC, halo_dir::count, remote_spec.data());
}

void tausch_exchanger::exchange_func(int k, real_t * gf)
{
}

void tausch_exchanger::exchange_sten(int k, real_t * so)
{
}

void tausch_exchanger::exchange(mpi::grid_func & f)
{
	std::array<int, halo_dir::count> tags;
	for (int dir = 0; dir < halo_dir::count; dir++)
		tags[dir] = dir;

	tausch->postAllReceives(TAUSCH_CwC, tags.data());

	// int rank;
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// printf("[%d] %d\n", rank, mpitag);

	for (int dir = 0; dir < halo_dir::count; dir++) {
		tausch->packSendBuffer(TAUSCH_CwC, dir, 0, f.data());
		tausch->send(TAUSCH_CwC, dir, tags[dir]);

		tausch->recv(TAUSCH_CwC, dir);
		tausch->unpackRecvBuffer(TAUSCH_CwC, dir, 0, f.data());
	}
}


}}}
