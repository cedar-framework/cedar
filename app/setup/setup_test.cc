#include <iostream>
#include <memory>

#include <boxmg-common.h>
#include <boxmg-2d.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;
	using namespace boxmg::bmg2d;

	len_t nx = 5;
	len_t ny = nx;

	auto so = StencilOp(nx,ny);
	GridStencil & sten = so.stencil();
	sten.five_pt() = true;

	for (auto j: sten.range(1)) {
		for (auto i: sten.range(0)) {
			sten(i,j,Dir::E) = 1;
			sten(i,j,Dir::N) = 1;
			sten(i,j,Dir::C) = 4;
			sten(i,j,Dir::S) = 1;
			sten(i,j,Dir::W) = 1;
		}
	}

	for (auto idx: sten.boarder(Dir::N)) sten(idx, Dir::N) = 0;
	for (auto idx: sten.boarder(Dir::S)) sten(idx, Dir::S) = 0;
	for (auto idx: sten.boarder(Dir::E)) sten(idx, Dir::E) = 0;
	for (auto idx: sten.boarder(Dir::W)) sten(idx, Dir::W) = 0;

	solver::BoxMG bmg(std::move(so));

	std::ofstream rfile;
	rfile.open("Restrict", std::ios::out | std::ios::trunc | std::ios::binary);
	rfile << bmg.level(bmg.nlevels()-1).P;
	rfile.close();

	std::ofstream sten_file;
	sten_file.open("Stencil", std::ios::out | std::ios::trunc | std::ios::binary);
	sten_file << bmg.level(bmg.nlevels()-1).A;
	sten_file.close();

	log::info << "Finished Test" << std::endl;

	return 0;
}
