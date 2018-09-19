#include <cedar/3d/util/io.h>

namespace cedar { namespace cdr3 { namespace util {

template<class sten>
void writeascii(mpi::stencil_op<sten> & so)
{
	auto & topo = so.grid();
	std::size_t nvals = (topo.nlocal(0)+1) * (topo.nlocal(1)+1) * (topo.nlocal(2)+1) * stencil_ndirs<sten>::value;

	std::string fname("so-" + std::to_string(topo.coord(0)) + "."
	                  + std::to_string(topo.coord(1)) + "."
	                  + std::to_string(topo.coord(2)));
	std::ofstream ofile(fname);

	ofile << stencil_ndirs<sten>::value << '\n';
	ofile << nvals << '\n';

	real_t * curr = so.data();
	for (auto i : range(nvals)) {
		(void)i;
		ofile << *curr << '\n';
		++curr;
	}
}


template void writeascii<xxvii_pt>(mpi::stencil_op<xxvii_pt> & so);
template void writeascii<seven_pt>(mpi::stencil_op<seven_pt> & so);

template<class sten>
void readascii(mpi::stencil_op<sten> & so)
{
	auto & topo = so.grid();
	std::size_t nvals = (topo.nlocal(0)+1) * (topo.nlocal(1)+1) * (topo.nlocal(2)+1) * stencil_ndirs<sten>::value;
	std::string fname("so-" + std::to_string(topo.coord(0)) + "."
	                  + std::to_string(topo.coord(1)) + "."
	                  + std::to_string(topo.coord(2)));
	std::ifstream ifile(fname);
	std::string line;

	std::size_t nvals_in;
	std::getline(ifile, line);
	std::getline(ifile, line);
	std::stringstream ss_nvals(line);
	ss_nvals >> nvals_in;

	if (nvals != nvals_in) {
		log::error << "input file nvals != stencil op nvals!" << std::endl;
	}

	real_t * curr = so.data();
	while (std::getline(ifile, line)) {
		std::stringstream ss(line);
		ss >> *curr;
		++curr;
	}
}


template void readascii<seven_pt>(mpi::stencil_op<seven_pt> & so);
template void readascii<xxvii_pt>(mpi::stencil_op<xxvii_pt> & so);

}}}
