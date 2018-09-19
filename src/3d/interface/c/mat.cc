#include <cedar/3d/mpi/stencil_op.h>

#include <cedar/capi.h>

#include <cedar/interface/object.h>
#include <cedar/interface/topo.h>
#include <cedar/interface/vec.h>
#include <cedar/interface/mat.h>


extern "C"
{
	int cedar_mat_create3d(cedar_topo topo, cedar_stencil_3d sten, cedar_mat *mat)
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		std::shared_ptr<grid_topo> topo_obj = cedar_topo_getobj(topo);
		if (not topo_obj) {
			*mat = CEDAR_MAT_NULL;
			return CEDAR_ERR_TOPO;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_MAT);

		cedar_mat_cont *cont = new cedar_mat_cont;
		cont->nd = 3;
		cont->topo = topo;
		cedar_object_incref(topo);
		config conf("config.json");
		cont->kman3 = mpi::build_kernel_manager(conf);
		cont->is_halo_setup = false;
		if (sten == CEDAR_STENCIL_SEVEN_PT) {
			cont->op3comp = std::make_unique<mpi::stencil_op<seven_pt>>(topo_obj);
			cont->compressed = true;
		} else {
			cont->op3full = std::make_unique<mpi::stencil_op<xxvii_pt>>(topo_obj);
			cont->compressed = false;
		}
		obj->ptr = reinterpret_cast<void*>(cont);
		*mat = obj->handle;

		return CEDAR_SUCCESS;
	}


	int cedar_mat_set3d(cedar_mat mat, unsigned int nvals, cedar_coord_3d coords[], cedar_real vals[])
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		auto *cont = cedar_mat_getobj(mat);
		if (not cont)
			return CEDAR_ERR_MAT;
		if (cont->nd != 3)
			return CEDAR_ERR_DIM;

		if (cont->compressed) {
			auto & so = *cont->op3comp;
			auto & grid = so.grid();
			for (auto i : range(nvals)) {
				len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
				len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
				len_t ck = static_cast<len_t>(coords[i].k - grid.is(2) + 2);
				// boxmg likes positive stencil coefficients
				if (coords[i].dir != 0) vals[i] = -1*vals[i];

				so(static_cast<len_t>(ci),
				   static_cast<len_t>(cj),
				   static_cast<len_t>(ck),
				   static_cast<seven_pt>(coords[i].dir)) = vals[i];
			}
		} else {
			auto & so = *cont->op3full;
			auto & grid = so.grid();
			for (auto i : range(nvals)) {
				len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
				len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
				len_t ck = static_cast<len_t>(coords[i].k - grid.is(2) + 2);
				// boxmg likes positive stencil coefficients
				if (coords[i].dir != 0) vals[i] = -1*vals[i];

				so(static_cast<len_t>(ci),
				   static_cast<len_t>(cj),
				   static_cast<len_t>(ck),
				   static_cast<xxvii_pt>(coords[i].dir)) = vals[i];
			}
		}

		return CEDAR_SUCCESS;
	}
}
