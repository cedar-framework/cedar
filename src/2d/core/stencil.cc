#include "stencil.h"

#include "types.h"

using namespace boxmg;
using namespace boxmg::bmg2d::core;

const real_t & Stencil::C() const { return vals[static_cast<int>(Dir::C)]; };
	real_t & Stencil::C() { return vals[static_cast<int>(Dir::C)]; };


	const real_t & Stencil::N() const { return vals[static_cast<int>(Dir::N)]; };
	real_t & Stencil::N() { return vals[static_cast<int>(Dir::N)]; };


	const real_t & Stencil::NE() const { return vals[static_cast<int>(Dir::NE)]; };
	real_t & Stencil::NE() { return vals[static_cast<int>(Dir::NE)]; };


	const real_t & Stencil::E() const { return vals[static_cast<int>(Dir::E)]; };
	real_t & Stencil::E() { return vals[static_cast<int>(Dir::E)]; };


	const real_t & Stencil::SE() const { return vals[static_cast<int>(Dir::SE)]; };
	real_t & Stencil::SE() { return vals[static_cast<int>(Dir::SE)]; };


	const real_t & Stencil::S() const { return vals[static_cast<int>(Dir::S)]; };
	real_t & Stencil::S() { return vals[static_cast<int>(Dir::S)]; };


	const real_t & Stencil::SW() const { return vals[static_cast<int>(Dir::SW)]; };
	real_t & Stencil::SW() { return vals[static_cast<int>(Dir::SW)]; };


	const real_t & Stencil::W() const { return vals[static_cast<int>(Dir::W)]; };
	real_t & Stencil::W() { return vals[static_cast<int>(Dir::W)]; };


	const real_t & Stencil::NW() const { return vals[static_cast<int>(Dir::NW)]; };
	real_t & Stencil::NW() { return vals[static_cast<int>(Dir::NW)]; };


namespace boxmg { namespace bmg2d { namespace core {
std::ostream & operator <<(std::ostream &os, const Stencil &sten)
{
	int width = 6;

		os << std::setw(width) << sten.NW() << " "
		   << std::setw(width) << sten.N()  << " "
		   << std::setw(width) << sten.NE() << std::endl;
		os << std::setw(width) << sten.W()  << " "
		   << std::setw(width) << sten.C()  << " "
		   << std::setw(width) << sten.E()  << std::endl;
		os << std::setw(width) << sten.SW() << " "
		   << std::setw(width) << sten.S()  << " "
		   << std::setw(width) << sten.SE() << std::endl;

		return os;
}
}}}
