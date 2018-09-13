#include <cedar/2d/mpi/solver.h>
#ifdef ENABLE_3D
#include <cedar/3d/mpi/solver.h>
#endif

#include <cedar/capi.h>
#include <cedar/interface/object.h>
#include <cedar/interface/vec.h>
#include <cedar/interface/mat.h>
#include <cedar/interface/solver.h>


cedar_solver_cont *cedar_solver_getobj(cedar_solver handle)
{
	if ((handle & CEDAR_KIND_MASK) != CEDAR_KIND_SOLVER)
		return nullptr;

	cedar_object *obj;
	int ierr = cedar_object_get(handle, &obj);
	if (ierr or (not obj))
		return nullptr;

	return reinterpret_cast<cedar_solver_cont*>(obj->ptr);
}


extern "C"
{
	int cedar_solver_create(cedar_mat mat, cedar_solver *solver)
	{
		using namespace cedar;

		auto *matobj = cedar_mat_getobj(mat);
		if (not matobj) {
			*solver = CEDAR_SOLVER_NULL;
			return CEDAR_ERR_MAT;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_SOLVER);

		cedar_solver_cont *slvobj = new cedar_solver_cont;
		slvobj->nd = matobj->nd;
		slvobj->compressed = matobj->compressed;

		if (matobj->nd == 2) {
			if (matobj->compressed)
				slvobj->slv2comp = std::make_unique<cdr2::mpi::solver<cdr2::five_pt>>(*(matobj->op2comp));
			else
				slvobj->slv2full = std::make_unique<cdr2::mpi::solver<cdr2::nine_pt>>(*(matobj->op2full));
		} else if (matobj->nd == 3) {
			#ifdef ENABLE_3D
			if (matobj->compressed)
				slvobj->slv3comp = std::make_unique<cdr3::mpi::solver<cdr3::seven_pt>>(*(matobj->op3comp));
			else
				slvobj->slv3full = std::make_unique<cdr3::mpi::solver<cdr3::xxvii_pt>>(*(matobj->op3full));
			#else
			return CEDAR_ERR_DIM;
			#endif
		}

		obj->ptr = reinterpret_cast<void*>(slvobj);
		*solver = obj->handle;

		return CEDAR_SUCCESS;
	}


	int cedar_solver_run(cedar_solver slv, cedar_vec x, cedar_vec b)
	{
		auto *slvobj = cedar_solver_getobj(slv);
		if (not slvobj)
			return CEDAR_ERR_SOLVER;

		auto *xobj = cedar_vec_getobj(x);
		if (not xobj)
			return CEDAR_ERR_VEC;

		auto *bobj = cedar_vec_getobj(b);
		if (not bobj)
			return CEDAR_ERR_VEC;

		if ((slvobj->nd != bobj->nd) or (slvobj->nd != xobj->nd))
			return CEDAR_ERR_DIM;

		if (slvobj->nd == 2) {
			if (slvobj->compressed)
				slvobj->slv2comp->solve(*(xobj->gfunc2), *(bobj->gfunc2));
			else
				slvobj->slv2full->solve(*(xobj->gfunc2), *(bobj->gfunc2));
		} else {
			#ifdef ENABLE_3D
			if (slvobj->compressed)
				slvobj->slv3comp->solve(*(xobj->gfunc3), *(bobj->gfunc3));
			else
				slvobj->slv3full->solve(*(xobj->gfunc3), *(bobj->gfunc3));
			#else
			return CEDAR_ERR_DIM;
			#endif
		}

		return CEDAR_SUCCESS;
	}


	int cedar_solver_free(cedar_solver *slv)
	{
		if ((*slv & CEDAR_KIND_MASK) != CEDAR_KIND_SOLVER)
			return CEDAR_ERR_SOLVER;

		cedar_object *obj;
		cedar_object_get(*slv, &obj);
		if (obj) {
			auto *cont = reinterpret_cast<cedar_solver_cont*>(obj->ptr);
			delete cont;
			cedar_object_free(*slv);
		} else
			return CEDAR_ERR_SOLVER;

		*slv = CEDAR_SOLVER_NULL;
		return CEDAR_SUCCESS;
	}
}
