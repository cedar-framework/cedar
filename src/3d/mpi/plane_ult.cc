#include <cedar/3d/mpi/plane_ult.h>

namespace cedar { namespace cdr3 { namespace mpi {

void plane_ult_run_full(void *args)
{
	ult_params<cdr2::nine_pt> *params = (ult_params<cdr2::nine_pt>*) args;
	auto & slv = *params->slv;

	auto & x = slv.levels.get(0).x;
	auto & b = slv.levels.get(0).b;

	slv.solve(b, x);
}


void plane_ult_run_comp(void *args)
{
	ult_params<cdr2::five_pt> *params = (ult_params<cdr2::five_pt>*) args;
	auto & slv = *params->slv;

	// auto & x = slv.levels.get(1).x;
	// std::size_t nbytes = x.len(0) * x.len(1) * sizeof(real_t);
	// real_t * xaddr = x.data();
	// printf("%p %p\n", xaddr, xaddr + (x.len(0) * x.len(1)));

	auto & x = slv.levels.template get<cdr2::five_pt>(0).x;
	auto & b = slv.levels.template get<cdr2::five_pt>(0).b;

	slv.solve(b, x);
}

}}}
