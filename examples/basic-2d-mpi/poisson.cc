#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <cedar/types.h>
#include <cedar/2d/mpi/grid_func.h>
#include <cedar/2d/util/topo.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/2d/mpi/gallery.h>

#include <cedar/util/time_log.h>

#include <ftl/Cedar.hpp>

static void set_problem(cedar::cdr2::mpi::grid_func & b)
{
	using namespace cedar;
	using namespace cedar::cdr2;

	const double pi = M_PI;

	auto rhs = [pi](real_t x, real_t y) {
		return 8*(pi*pi)*sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = b.grid();

	b.set(0);

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	real_t h2 = hx*hy;

	real_t nlx = topo.nlocal(0) - 2;
	real_t nly = topo.nlocal(1) - 2;

	real_t i1 = nlx + 1;
	real_t j1 = nly + 1;

	for (auto j : range<len_t>(1, j1)) {
		for (auto i : range<len_t>(1, i1)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			b(i,j) = rhs(x, y) * h2;

		}
	}
}


static void set_solution(cedar::cdr2::mpi::grid_func & q)
{
	using namespace cedar;

	const double pi = M_PI;

	auto sol = [pi](real_t x, real_t y) {
		return sin(2*pi*x)*sin(2*pi*y);
	};

	auto & topo = q.grid();

	real_t igs = topo.is(0);
	real_t jgs = topo.is(1);

	real_t hx = 1.0 / (topo.nglobal(0) - 1);
	real_t hy = 1.0 / (topo.nglobal(1) - 1);

	for (auto j : q.range(1)) {
		for (auto i : q.range(0)) {
			len_t is = igs + i;
			len_t js = jgs + j;

			real_t x = (is-1)*hx;
			real_t y = (js-1)*hy;

			q(i,j) = sol(x,y);
		}
	}
}


template <typename dtype, unsigned short ND>
ftl::Buffer<dtype> to_buffer(cedar::array<dtype, ND>& array) {
    std::vector<int32_t> shape(ND);

    for (std::size_t i = 0; i < ND; ++i) {
        const int32_t dim = array.len(i);
        shape[i] = dim;
    }

    ftl::Buffer<dtype> buf(array.data(), shape);

    return buf;
}


int main(int argc, char *argv[])
{
	using namespace cedar;
	using namespace cedar::cdr2;

	int provided;

	MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	timer_init(MPI_COMM_WORLD);

	config conf;

	auto grid = util::create_topo(conf);

	auto so = mpi::gallery::poisson(grid);
	mpi::grid_func b(grid);

        // auto* scp = so.data();
        // for (std::size_t i = 0; i < 7 * 7 * 3; ++i) {
        //     std::cerr << scp[i] << " ";
        // }
        // std::cerr << std::endl;

	set_problem(b);

	mpi::solver<five_pt> bmg(so);
        const std::size_t levels = bmg.nlevels();
        std::cout << "Solver has " << levels << " levels." << std::endl;

        for (std::size_t i = 0; i < levels; ++i) {
            if (i == 0) {
                auto level = bmg.levels.get<five_pt>(i);
                auto so = level.A;
                std::cerr << i << " " << so.len(0) << " " << so.len(1) << std::endl;
            } else {
                auto level = bmg.levels.get<nine_pt>(i);
                auto so = level.A;
                std::cerr << i << " " << so.len(0) << " " << so.len(1) << std::endl;
            }
        }

        // auto ci_buf = to_buffer(ci);
        // std::cout << ci_buf << std::endl;

        MPI_Barrier(MPI_COMM_WORLD);

        auto level = bmg.levels.get<five_pt>(0);
        auto ci = level.P;
        std::cout << ci.len(0) << " " << ci.len(1) << " " << ci.len(2) << std::endl;


	MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
	auto sol = bmg.solve(b);

        // auto level = bmg.levels.get<five_pt>(0);
        // auto ci = level.P;
        ci = level.P;
        std::cout << ci.len(0) << " " << ci.len(1) << " " << ci.len(2) << std::endl;


	// mpi::grid_func exact_sol(sol.grid_ptr());
	// set_solution(exact_sol);

	// mpi::grid_func diff = exact_sol - sol;

	//log::status << "Solution norm: " << diff.inf_norm() << std::endl;

	//timer_save("timings.json");

	//log::status << "Finished Test" << std::endl;

        if (rank == size - 1) {
            std::cout << std::endl << "Press enter to exit..." << std::endl;
            std::cin.get();
        } else {
            std::cout << std::endl << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
