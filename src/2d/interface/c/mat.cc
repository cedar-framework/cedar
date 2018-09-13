#include <cedar/2d/mpi/stencil_op.h>

#include <cedar/capi.h>

#include <cedar/interface/object.h>
#include <cedar/interface/topo.h>
#include <cedar/interface/vec.h>
#include <cedar/interface/mat.h>


cedar_mat_cont *cedar_mat_getobj(cedar_mat handle)
{
	if ((handle & CEDAR_KIND_MASK) != CEDAR_KIND_MAT)
		return nullptr;

	cedar_object *obj;
	int ierr = cedar_object_get(handle, &obj);
	if (ierr or (not obj))
		return nullptr;

	return reinterpret_cast<cedar_mat_cont*>(obj->ptr);
}


extern "C"
{
	int cedar_mat_create2d(cedar_topo topo, cedar_stencil_2d sten, cedar_mat *mat)
	{
		using namespace cedar;
		using namespace cedar::cdr2;

		std::shared_ptr<grid_topo> topo_obj = cedar_topo_getobj(topo);
		if (not topo_obj) {
			*mat = CEDAR_MAT_NULL;
			return CEDAR_ERR_TOPO;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_MAT);

		cedar_mat_cont *cont = new cedar_mat_cont;
		cont->nd = 2;
		cont->topo = topo;
		cedar_object_incref(topo);
		if (sten == CEDAR_STENCIL_FIVE_PT) {
			cont->op2comp = std::make_unique<mpi::stencil_op<five_pt>>(topo_obj);
			cont->compressed = true;
		} else {
			cont->op2full = std::make_unique<mpi::stencil_op<nine_pt>>(topo_obj);
			cont->compressed = false;
		}
		obj->ptr = reinterpret_cast<void*>(cont);
		*mat = obj->handle;

		return CEDAR_SUCCESS;
	}


	int cedar_mat_set2d(cedar_mat mat, unsigned int nvals, cedar_coord_2d coords[], cedar_real vals[])
	{
		using namespace cedar;
		using namespace cedar::cdr2;

		auto *cont = cedar_mat_getobj(mat);
		if (not cont)
			return CEDAR_ERR_MAT;
		if (cont->nd != 2)
			return CEDAR_ERR_DIM;

		if (cont->compressed) {
			auto & so = *cont->op2comp;
			auto & grid = so.grid();
			for (auto i : range(nvals)) {
				len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
				len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
				auto dir = coords[i].dir;
				// boxmg likes positive stencil coefficients
				if (dir != 0) vals[i] = -1*vals[i];

				// vertex based input -> boxmg symmetric storage
				if (dir == BMG2_SE) {
					return CEDAR_ERR_STENCIL;
				} else if (dir == BMG2_N) {
					cj++;
					dir = BMG2_S;
				} else if (dir == BMG2_NE) {
					return CEDAR_ERR_STENCIL;
				} else if (dir == BMG2_E) {
					ci++;
					dir = BMG2_W;
				} else if (dir == BMG2_NW)
					return CEDAR_ERR_STENCIL;

				so(ci, cj, static_cast<five_pt>(dir)) = vals[i];
			}
		} else {
			auto & so = *cont->op2full;
			auto & grid = so.grid();
			for (auto i : range(nvals)) {
				len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
				len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
				auto dir = coords[i].dir;
				// boxmg likes positive stencil coefficients
				if (dir != 0) vals[i] = -1*vals[i];

				// vertex based input -> boxmg symmetric storage
				if (dir == BMG2_SE) {
					ci++;
					dir = BMG2_NW;
				} else if (dir == BMG2_N) {
					cj++;
					dir = BMG2_S;
				} else if (dir == BMG2_NE) {
					ci++; cj++;
					dir = BMG2_SW;
				} else if (dir == BMG2_E) {
					ci++;
					dir = BMG2_W;
				} else if (dir == BMG2_NW)
					cj++;

				so(ci, cj, static_cast<nine_pt>(dir)) = vals[i];
			}
		}

		return CEDAR_SUCCESS;
	}


	int cedar_mat_free(cedar_mat *mat)
	{
		if ((*mat & CEDAR_KIND_MASK) != CEDAR_KIND_MAT)
			return CEDAR_ERR_MAT;

		cedar_object *obj;
		cedar_object_get(*mat, &obj);
		if (obj) {
			auto *cont = reinterpret_cast<cedar_mat_cont*>(obj->ptr);
			cedar_object_decref(cont->topo);
			delete cont;
			cedar_object_free(*mat);
		} else
			return CEDAR_ERR_MAT;

		*mat = CEDAR_MAT_NULL;
		return CEDAR_SUCCESS;
	}
}
