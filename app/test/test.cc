#include <mpi.h>
#include <iostream>
#include <memory>
#include <math.h>

#include <boxmg/perf/params.h>
#include <boxmg/types.h>


int main(int argc, char *argv[])
{
	using namespace boxmg;

	config::reader conf;

	log::status << params::compute_tc(2, conf) << std::endl;;

	return 0;
}
