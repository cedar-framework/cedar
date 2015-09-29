#include "kernel/manager.h"
#include "kernel/name.h"

#include "registry.h"

using namespace boxmg::bmg2d::kernel;

void Registry::setup_interp(int kf, int kc, int nog, const core::StencilOp & fop,
                            const core::StencilOp &cop, inter::ProlongOp & P)
{
	active.run(name::setup_interp,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             fop,cop,P);
}


void Registry::galerkin_prod(int kf, int kc, int nog,
                   const inter::ProlongOp & P,
                   const core::StencilOp & fop,
                   core::StencilOp & cop)
{
	active.run(name::galerkin_prod,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             P,fop,cop);


}


void Registry::setup_relax(const core::StencilOp & so,
                 core::RelaxStencil & sor)
{
	active.run(name::setup_relax, so, sor);
}


void Registry::setup_relax_x(const core::StencilOp & so,
                             core::RelaxStencil & sor)
{
	active.run(name::setup_relax_x, so, sor);
}


void Registry::setup_relax_y(const core::StencilOp & so,
                             core::RelaxStencil & sor)
{
	active.run(name::setup_relax_y, so, sor);
}


void Registry::setup_cg_lu(const core::StencilOp & so,
                 core::GridFunc & abd)
{
	active.run(name::setup_cg_lu, so, abd);
}

void Registry::relax(const core::StencilOp & so,
           core::GridFunc & x,
           const core::GridFunc & b,
           const core::RelaxStencil & sor,
           cycle::Dir cycle_dir)
{
	active.run(name::relax, so, x, b, sor, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::relax_lines_x(const core::StencilOp & so,
           core::GridFunc & x,
           const core::GridFunc & b,
           const core::RelaxStencil & sor,
           cycle::Dir cycle_dir)
{
	active.run(name::relax_lines_x, so, x, b, sor, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::relax_lines_y(const core::StencilOp & so,
           core::GridFunc & x,
           const core::GridFunc & b,
           const core::RelaxStencil & sor,
           cycle::Dir cycle_dir)
{
	active.run(name::relax_lines_y, so, x, b, sor, static_cast<cycle::Dir>(cycle_dir));
}


void Registry::solve_cg(core::GridFunc &x,
                        const core::GridFunc &b,
                        const core::GridFunc &ABD,
                        real_t *bbd)
{
	active.run(name::solve_cg, x, b, ABD, static_cast<real_t*>(bbd));
}

void Registry::setup_nog(core::mpi::GridTopo &topo,
               len_t min_coarse, int *nog)
{
	active.run(name::setup_nog, topo,
	             static_cast<len_t>(min_coarse),
	             static_cast<int*>(nog));
}


void Registry::halo_setup(core::mpi::GridTopo &topo,
                void **halo_ctx)
{
	active.run(name::halo_setup,
	             static_cast<core::mpi::GridTopo&>(topo),
	             static_cast<void**>(halo_ctx));
}

void Registry::halo_exchange(core::mpi::GridFunc & f)
{
	active.run(name::halo_exchange,
	             static_cast<core::mpi::GridFunc&>(f));
}


void Registry::halo_exchange(const core::mpi::GridFunc & f, void *halo_ctx)
{
	auto & fd = const_cast<core::mpi::GridFunc&>(f);
	fd.halo_ctx = halo_ctx;
	active.run(name::halo_exchange,
	           static_cast<core::mpi::GridFunc&>(fd));
}


void Registry::halo_stencil_exchange(core::mpi::StencilOp & so)
{
	active.run(name::halo_stencil_exchange,
	             so);
}


void Registry::setup_cg_boxmg(const core::StencilOp & so,
                              std::shared_ptr<solver::BoxMG> *solver)
{
	active.run(name::setup_cg_boxmg, so,
	           static_cast<std::shared_ptr<solver::BoxMG>*>(solver));
}


void Registry::solve_cg_boxmg(const solver::BoxMG &cg_solver,
                              core::GridFunc &x,
                              const core::GridFunc &b)
{
	active.run(name::solve_cg_boxmg, cg_solver, x, b);
}


void Registry::matvec(const core::StencilOp & so,
                      const core::GridFunc & x,
                      core::GridFunc &b)
{
	active.run(name::matvec, so, x, b);
}
