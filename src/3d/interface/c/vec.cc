#include <cedar/3d/mpi/grid_func.h>

#include <cedar/capi.h>

#include <cedar/interface/object.h>
#include <cedar/interface/topo.h>
#include <cedar/interface/vec.h>


extern "C"
{
	int cedar_vec_create3d(cedar_topo topo, cedar_vec *vec)
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		std::shared_ptr<grid_topo> topo_obj = cedar_topo_getobj(topo);
		if (not topo_obj) {
			*vec = CEDAR_VEC_NULL;
			return CEDAR_ERR_TOPO;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_VEC);

		cedar_vec_cont *cont = new cedar_vec_cont;
		cont->nd = 3;
		cont->topo = topo;
		cedar_object_incref(topo);
		cont->gfunc3 = std::make_unique<mpi::grid_func>(topo_obj);
		obj->ptr = reinterpret_cast<void*>(cont);
		*vec = obj->handle;

		return CEDAR_SUCCESS;
	}
}
