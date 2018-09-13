#include <cedar/2d/mpi/grid_func.h>

#include <cedar/capi.h>

#include <cedar/interface/object.h>
#include <cedar/interface/topo.h>
#include <cedar/interface/vec.h>


cedar_vec_cont *cedar_vec_getobj(cedar_vec handle)
{
	if ((handle & CEDAR_KIND_MASK) != CEDAR_KIND_VEC)
		return nullptr;

	cedar_object *obj;
	int ierr = cedar_object_get(handle, &obj);
	if (ierr or (not obj))
		return nullptr;

	return reinterpret_cast<cedar_vec_cont*>(obj->ptr);
}


extern "C"
{
	int cedar_vec_create2d(cedar_topo topo, cedar_vec *vec)
	{
		using namespace cedar;
		using namespace cedar::cdr2;

		std::shared_ptr<grid_topo> topo_obj = cedar_topo_getobj(topo);
		if (not topo_obj) {
			*vec = CEDAR_VEC_NULL;
			return CEDAR_ERR_TOPO;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_VEC);

		cedar_vec_cont *cont = new cedar_vec_cont;
		cont->nd = 2;
		cont->topo = topo;
		cedar_object_incref(topo);
		cont->gfunc2 = std::make_unique<mpi::grid_func>(topo_obj);
		obj->ptr = reinterpret_cast<void*>(cont);
		*vec = obj->handle;

		return CEDAR_SUCCESS;
	}


	int cedar_vec_free(cedar_vec *vec)
	{
		if ((*vec & CEDAR_KIND_MASK) != CEDAR_KIND_TOPO)
			return CEDAR_ERR_VEC;

		cedar_object *obj;
		cedar_object_get(*vec, &obj);
		if (obj) {
			auto *cont = reinterpret_cast<cedar_vec_cont*>(obj->ptr);
			cedar_object_decref(cont->topo);
			delete cont;
			cedar_object_free(*vec);
		} else
			return CEDAR_ERR_VEC;

		*vec = CEDAR_VEC_NULL;
		return CEDAR_SUCCESS;
	}
}
