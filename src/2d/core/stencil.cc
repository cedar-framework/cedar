#include "stencil.h"

#include "types.h"

using namespace boxmg;
using namespace boxmg::bmg2d::core;

const real_t & Stencil::C() const { return vals[static_cast<int>(dir::C)]; };
	real_t & Stencil::C() { return vals[static_cast<int>(dir::C)]; };


	const real_t & Stencil::N() const { return vals[static_cast<int>(dir::N)]; };
	real_t & Stencil::N() { return vals[static_cast<int>(dir::N)]; };


	const real_t & Stencil::NE() const { return vals[static_cast<int>(dir::NE)]; };
	real_t & Stencil::NE() { return vals[static_cast<int>(dir::NE)]; };


	const real_t & Stencil::E() const { return vals[static_cast<int>(dir::E)]; };
	real_t & Stencil::E() { return vals[static_cast<int>(dir::E)]; };


	const real_t & Stencil::SE() const { return vals[static_cast<int>(dir::SE)]; };
	real_t & Stencil::SE() { return vals[static_cast<int>(dir::SE)]; };


	const real_t & Stencil::S() const { return vals[static_cast<int>(dir::S)]; };
	real_t & Stencil::S() { return vals[static_cast<int>(dir::S)]; };


	const real_t & Stencil::SW() const { return vals[static_cast<int>(dir::SW)]; };
	real_t & Stencil::SW() { return vals[static_cast<int>(dir::SW)]; };


	const real_t & Stencil::W() const { return vals[static_cast<int>(dir::W)]; };
	real_t & Stencil::W() { return vals[static_cast<int>(dir::W)]; };


	const real_t & Stencil::NW() const { return vals[static_cast<int>(dir::NW)]; };
	real_t & Stencil::NW() { return vals[static_cast<int>(dir::NW)]; };


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
