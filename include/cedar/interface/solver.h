#ifndef CEDAR_INTERFACE_SOLVER_H
#define CEDAR_INTERFACE_SOLVER_H

#include <memory>

#include <cedar/capi.h>
#include <cedar/2d/mpi/solver.h>
#ifdef ENABLE_3D
#include <cedar/3d/mpi/solver.h>
#endif

struct cedar_solver_cont
{
	unsigned short nd;
	bool compressed;
	std::unique_ptr<cedar::cdr2::mpi::solver<cedar::cdr2::five_pt>> slv2comp;
	std::unique_ptr<cedar::cdr2::mpi::solver<cedar::cdr2::nine_pt>> slv2full;
	#ifdef ENABLE_3D
	std::unique_ptr<cedar::cdr3::mpi::solver<cedar::cdr3::seven_pt>> slv3comp;
	std::unique_ptr<cedar::cdr3::mpi::solver<cedar::cdr3::xxvii_pt>> slv3full;
	#endif
};


cedar_solver_cont *cedar_solver_getobj(cedar_solver slv);

#endif
