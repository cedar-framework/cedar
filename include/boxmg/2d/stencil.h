#ifndef BOXMG_2D_CORE_STENCIL_H
#define BOXMG_2D_CORE_STENCIL_H

#include <iostream>
#include <iomanip>


namespace boxmg {
namespace bmg2d {
namespace core {

    class Stencil
    {
	public:

		Stencil() {};

		friend std::ostream &operator << (std::ostream &os, const Stencil &sten);

		const real_t & C() const;
		real_t & C();

		const real_t & N() const;
		real_t & N();

		const real_t & NE() const;
		real_t & NE();

		const real_t & E() const;
		real_t & E();

		const real_t & SE() const;
		real_t & SE();

		const real_t & S() const;
		real_t & S();

		const real_t & SW() const;
		real_t & SW();

		const real_t & W() const;
		real_t & W();

		const real_t & NW() const;
		real_t & NW();


	private:
		real_t vals[9];
    };
}
}
}

#endif
