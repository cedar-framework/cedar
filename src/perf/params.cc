#include <mpi.h>
#include <random>
#include <boxmg/2d/grid_func.h>
#include <boxmg/2d/stencil_op.h>
#include <boxmg/2d/kernel/factory.h>
#include <boxmg/perf/params.h>

using namespace boxmg;


static void fill_random(real_t * arr, len_t size)
{
	for (len_t i = 0; i < size; i++) {
		arr[i] = static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX);
	}
}


float params::compute_tc(int nd, config::reader & conf)
{
	// float tc;
	// int nsamples = 100;
	// int nlocal = 40;

	// if (nd == 2) {
	// 	using namespace boxmg::bmg2d;
	// 	auto kreg = boxmg::bmg2d::kernel::factory::from_config(conf);
	// 	int ns = 9;
	// 	auto so = stencil_op(nlocal, nlocal);
	// 	so.set_registry(kreg);
	// 	auto & sten = so.stencil();
	// 	auto len = sten.len(0)*sten.len(1);
	// 	fill_random(sten.data(), len*5);
	// 	grid_func b(nlocal, nlocal);
	// 	grid_func x(nlocal, nlocal);
	// 	grid_func r(nlocal, nlocal);
	// 	fill_random(b.data(), len);
	// 	fill_random(x.data(), len);
	// 	double starttime, endtime;
	// 	starttime = MPI_Wtime();
	// 	for (int i = 0; i < nsamples; i++) {
	// 		so.residual(x, b, r);
	// 	}
	// 	endtime = MPI_Wtime();
	// 	tc = (endtime-starttime) / nsamples / 2 / ns / (sten.shape(0)*sten.shape(1));
	// } else {
	// 	int ns = 27;
	// }

	// return tc;
	return conf.get<float>("machine.fp_perf");
}
