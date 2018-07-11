#include <cedar/3d/mpi/plane_ult.h>

namespace cedar { namespace cdr3 { namespace mpi {

void plane_ult_run_full(void *args)
{
	ult_params<cdr2::nine_pt> *params = (ult_params<cdr2::nine_pt>*) args;
	auto & slv = *params->slv;
}


void plane_ult_run_comp(void *args)
{
	ult_params<cdr2::five_pt> *params = (ult_params<cdr2::five_pt>*) args;
	auto & slv = *params->slv;
}

}}}
