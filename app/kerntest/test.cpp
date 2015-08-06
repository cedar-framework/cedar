#include <iostream>
#include <memory>

#include <boxmg-common.h>
#include <boxmg-2d.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;
	using namespace boxmg::bmg2d::core;

	StencilOp so(4,5);
	GridStencil & sten = so.stencil();

	for (auto j: sten.range(0)) {
		for (auto i: sten.range(1)) {
			sten(j,i,Dir::C) = 8;
			sten(j,i,Dir::S) = 1;
			sten(j,i,Dir::W) = 1;
			sten(j,i,Dir::SE) = 1;
			sten(j,i,Dir::SW) = 1;
			sten(j,i,Dir::N) = 1;
			sten(j,i,Dir::E) = 1;
			sten(j,i,Dir::NE) = 1;
			sten(j,i,Dir::NW) = 1;
		}
	}

	GridFunc x = GridFunc::ones(4,5);
	GridFunc b = GridFunc::ones(4,5);
	GridFunc r = GridFunc::ones(4,5);

	x.scale(3.5);
	b.scale(1.2);
	r.scale(5);

	using resk_t = Kernel<const StencilOp&, const GridFunc&, const GridFunc &, GridFunc&>;
	kernel::Manager::reg(kernel::name::residual, "no-op",
						 resk_t([](const StencilOp&,
								   const GridFunc&,
								   const GridFunc&,
								   GridFunc&) -> void {
									std::cout << "No-Op" << std::endl;
								}));

	so.residual(x,b,r);

	std::cout << '\n' << r << std::endl;

	return 0;
}
