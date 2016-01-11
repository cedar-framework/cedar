#include <boxmg/types.h>
#include <boxmg/3d/grid_func.h>
#include <boxmg/3d/grid_stencil.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg3;

	grid_func v = grid_func::ones(2,2,2);
	v(2,1,2) = 3.0;

	log::status << "Writing Vector\n" << std::endl;
	std::cout << v << std::endl;
	log::status << v.lp_norm<2>() << std::endl;

	grid_stencil gs(3,3,3);
	gs(1,0,1,dir::P) = 3;
	log::status << gs(1,0,1,dir::P) << std::endl;

	return 0;
}
