#include "boxmg/2d/kernel/name.h"

#include "boxmg/2d/kernel/registry.h"

using namespace boxmg::bmg2d::kernel;

void Registry::setup_interp(int kf, int kc, int nog, const stencil_op & fop,
                            const stencil_op &cop, inter::prolong_op & P)
{
	active.run(name::setup_interp,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             fop,cop,P);
}


void Registry::galerkin_prod(int kf, int kc, int nog,
                   const inter::prolong_op & P,
                   const stencil_op & fop,
                   stencil_op & cop)
{
	active.run(name::galerkin_prod,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             P,fop,cop);


}


void Registry::setup_relax(const stencil_op & so,
                 relax_stencil & sor)
{
	active.run(name::setup_relax, so, sor);
}


void Registry::setup_relax_x(const stencil_op & so,
                             relax_stencil & sor)
{
	active.run(name::setup_relax_x, so, sor);
}


void Registry::setup_relax_y(const stencil_op & so,
                             relax_stencil & sor)
{
	active.run(name::setup_relax_y, so, sor);
}


void Registry::setup_cg_lu(const stencil_op & so,
                 grid_func & abd)
{
	active.run(name::setup_cg_lu, so, abd);
}

void Registry::relax(const stencil_op & so,
           grid_func & x,
           const grid_func & b,
           const relax_stencil & sor,
           cycle::Dir cycle_dir)
{
	active.run(name::relax, so, x, b, sor, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::relax_lines_x(const stencil_op & so,
           grid_func & x,
           const grid_func & b,
           const relax_stencil & sor,
           grid_func &res,
           cycle::Dir cycle_dir)
{
	active.run(name::relax_lines_x, so, x, b, sor, res, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::relax_lines_y(const stencil_op & so,
           grid_func & x,
           const grid_func & b,
           const relax_stencil & sor,
           grid_func &res,
           cycle::Dir cycle_dir)
{
	active.run(name::relax_lines_y, so, x, b, sor, res, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::solve_cg(grid_func &x,
                        const grid_func &b,
                        const grid_func &ABD,
                        real_t *bbd)
{
	active.run(name::solve_cg, x, b, ABD, static_cast<real_t*>(bbd));
}

void Registry::setup_nog(mpi::grid_topo &topo,
               len_t min_coarse, int *nog)
{
	active.run(name::setup_nog, topo,
	             static_cast<len_t>(min_coarse),
	             static_cast<int*>(nog));
}


void Registry::halo_setup(mpi::grid_topo &topo,
                void **halo_ctx)
{
	active.run(name::halo_setup,
	             static_cast<mpi::grid_topo&>(topo),
	             static_cast<void**>(halo_ctx));
}

void Registry::halo_exchange(mpi::grid_func & f)
{
	active.run(name::halo_exchange,
	             static_cast<mpi::grid_func&>(f));
}


void Registry::halo_exchange(const mpi::grid_func & f, void *halo_ctx)
{
	auto & fd = const_cast<mpi::grid_func&>(f);
	fd.halo_ctx = halo_ctx;
	active.run(name::halo_exchange,
	           static_cast<mpi::grid_func&>(fd));
}


void Registry::halo_stencil_exchange(mpi::stencil_op & so)
{
	active.run(name::halo_stencil_exchange,
	             so);
}


void Registry::setup_cg_boxmg(const stencil_op & so,
                              std::shared_ptr<solver> *bmg)
{
	active.run(name::setup_cg_boxmg, so,
	           static_cast<std::shared_ptr<solver>*>(bmg));
}


void Registry::solve_cg_boxmg(const solver &cg_solver,
                              grid_func &x,
                              const grid_func &b)
{
	active.run(name::solve_cg_boxmg, cg_solver, x, b);
}


void Registry::matvec(const stencil_op & so,
                      const grid_func & x,
                      grid_func &b)
{
	active.run(name::matvec, so, x, b);
}
