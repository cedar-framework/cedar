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
	}

	remote_spec[halo_dir::right].haloX = 0;
	remote_spec[halo_dir::right].haloY = 0;
	remote_spec[halo_dir::right].haloWidth = 1;
	remote_spec[halo_dir::right].haloHeight = topo.nlocal(1);
	remote_spec[halo_dir::right].remoteMpiRank = rank - 1;

	local_spec[halo_dir::right].haloX = topo.nlocal(0);
	local_spec[halo_dir::right].haloY = 0;
	local_spec[halo_dir::right].haloWidth = 1;
	local_spec[halo_dir::right].haloHeight = topo.nlocal(1);
	local_spec[halo_dir::right].remoteMpiRank = rank + 1;

	int nbuf = 1;

	tausch = std::make_unique<Tausch2D<double>>(MPI_DOUBLE, nbuf, nullptr, topo.comm);

	tausch->setLocalHaloInfo(TAUSCH_CwC, 1, local_spec.data());
	tausch->setRemoteHaloInfo(TAUSCH_CwC, 1, remote_spec.data());
}

void tausch_exchanger::exchange_func(int k, real_t * gf)
{
}

void tausch_exchanger::exchange_sten(int k, real_t * so)
{
}

void tausch_exchanger::exchange(mpi::grid_func & f)
{
	int mpitag = 0;
	tausch->postAllReceives(TAUSCH_CwC, &mpitag);

	tausch->packSendBuffer(TAUSCH_CwC, halo_dir::right, 0, f.data());
	tausch->send(TAUSCH_CwC, halo_dir::right, mpitag);

	tausch->recv(TAUSCH_CwC, halo_dir::right);
	tausch->unpackRecvBuffer(TAUSCH_CwC, halo_dir::right, 0, f.data());
}


}}}
