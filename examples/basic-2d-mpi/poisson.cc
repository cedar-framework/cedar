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
#include <ftl/Runtime.hpp>

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

    ftl::runtime_parse_args(argc, argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    timer_init(MPI_COMM_WORLD);

    config conf;
    std::shared_ptr<config> conf_ptr(&conf, [](config*){}); /* :-) */
    std::string timing_fname("timings.json");

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--n")) {
            int n = atoi(argv[++i]);
            conf.setvec<int>("grid.n", {n, n});
        } else if (!strcmp(argv[i], "--gpu")) {
            conf.set<bool>("solver.gpu", true);
        } else if (!strcmp(argv[i], "--cpu")) {
            conf.set<bool>("solver.gpu", false);
        } else if (!strcmp(argv[i], "--local")) {
            conf.set<bool>("grid.local", true);
        } else if (!strcmp(argv[i], "--global")) {
            conf.set<bool>("grid.local", false);
        } else if (!strcmp(argv[i], "--timing-out")) {
            timing_fname = std::string(argv[++i]);
        }
    }

    auto grid = util::create_topo(conf);
    auto so = mpi::gallery::poisson(grid);
    mpi::grid_func b(grid);

    std::cerr << "Running on " << b.len(0) << "x" << b.len(1) << " problem." << std::endl;

    set_problem(b);

    mpi::solver<five_pt> bmg(so, conf_ptr);
    const std::size_t levels = bmg.nlevels();

    MPI_Barrier(MPI_COMM_WORLD);

    auto level = bmg.levels.get<five_pt>(0);
    auto ci = level.P;

    MPI_Barrier(MPI_COMM_WORLD); // synchronize before timing solve
    auto sol = bmg.solve(b);

    timer_save(timing_fname);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
