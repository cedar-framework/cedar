#ifndef CEDAR_HALO_EXCHANGE_H
#define CEDAR_HALO_EXCHANGE_H

#include <cedar/mpi/grid_topo.h>
#include <cedar/array.h>
#include <cedar/service.h>

namespace cedar { namespace services {

class halo_exchange_base : public service
{
public:
	virtual void exchange_func(int k, real_t *gf) = 0;
	virtual void exchange_sten(int k, real_t *so) = 0;
};


template<class solver_types>
class halo_exchange : public halo_exchange_base
{
public:
	template<class sten>
	using stencil_op = typename solver_types::template stencil_op<sten>;
	using comp_sten = typename solver_types::comp_sten;
	using full_sten = typename solver_types::full_sten;
	using grid_func = typename solver_types::grid_func;

	const static std::string name() { return "halo exchange"; }

	virtual void setup(std::vector<topo_ptr> topos) = 0;
	virtual void run(stencil_op<comp_sten> & so) = 0;
	virtual void run(stencil_op<full_sten> & so) = 0;
	virtual void run(grid_func & gf, unsigned short dmask=7) = 0;

	virtual aarray<int, len_t, 2> & leveldims(int k) = 0;
	virtual len_t * datadist(int k, int grid) = 0;
	virtual std::vector<real_t> & linebuf() = 0;
};

}}

#endif
