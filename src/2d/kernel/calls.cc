#include "kernel/manager.h"
#include "kernel/name.h"
#include "kernel/calls.h"


namespace boxmg { namespace bmg2d { namespace kernel {

void setup_interp(int kf, int kc, int nog, const core::StencilOp & fop,
                  const core::StencilOp &cop, inter::ProlongOp & P)
{
	// Manager::run(name::setup_interp,kf,kc,nog,fop,cop,P);
	// Manager::run(name::setup_interp,
	//              static_cast<int>(kf),
	//              static_cast<int>(kc),
	//              static_cast<int>(nog),
	//              static_cast<const core::StencilOp&>(fop),
	//              static_cast<const core::StencilOp&>(cop),
	//              static_cast<inter::ProlongOp&>(P));
	Manager::run(name::setup_interp,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             fop,cop,P);
}


void galerkin_prod(int kf, int kc, int nog,
                   const inter::ProlongOp & P,
                   const core::StencilOp & fop,
                   core::StencilOp & cop)
{
	Manager::run(name::galerkin_prod,
	             static_cast<int>(kf),
	             static_cast<int>(kc),
	             static_cast<int>(nog),
	             P,fop,cop);


}


void setup_relax(const core::StencilOp & so,
                 core::RelaxStencil & sor)
{
	Manager::run(name::setup_relax, so, sor);
}


void setup_cg_lu(const core::StencilOp & so,
                 core::GridFunc & abd)
{
	Manager::run(name::setup_cg_lu, so, abd);
}

void relax(const core::StencilOp & so,
           core::GridFunc & x,
           const core::GridFunc & b,
           const core::RelaxStencil & sor,
           cycle::Dir cycle_dir)
{
	Manager::run(name::relax, so, x, b, sor, static_cast<cycle::Dir>(cycle_dir));
}


void solve_cg(core::GridFunc &x,
              const core::GridFunc &b,
              const core::GridFunc &ABD)
{
	Manager::run(name::solve_cg, x, b, ABD);
}

void setup_nog(core::mpi::GridTopo &topo,
               len_t min_coarse, int *nog)
{
	Manager::run(name::setup_nog, topo,
	             static_cast<len_t>(min_coarse),
	             static_cast<int*>(nog));
}


void halo_setup(core::mpi::GridTopo &topo,
                void **halo_ctx)
{
	Manager::run(name::halo_setup,
	             static_cast<core::mpi::GridTopo&>(topo),
	             static_cast<void**>(halo_ctx));
}

void halo_exchange(core::mpi::GridFunc & f)
{
	Manager::run(name::halo_exchange,
	             static_cast<core::mpi::GridFunc&>(f));
}


void halo_stencil_exchange(core::mpi::StencilOp & so)
{
	Manager::run(name::halo_stencil_exchange,
	             so);
}

}}}
