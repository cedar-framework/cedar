#include <stdbool.h>

#include <cedar/types.h>
#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/config/reader.h>

#include <cedar/3d/interface/c/operator.h>
#include <cedar/3d/interface/c/types.h>

extern "C"
{
	bmg3_operator bmg3_operator_create(bmg3_topo topo)
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		auto grid = *(reinterpret_cast<std::shared_ptr<cedar::grid_topo>*>(topo));

		config::reader conf("config.json");
		auto * op_cont = new op_container(grid, conf);

		return reinterpret_cast<bmg3_operator>(op_cont);
	}


	void bmg3_operator_set(bmg3_operator op, unsigned int nvals, grid_coord_3d coords[], double vals[])
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		auto * op_cont = reinterpret_cast<op_container*>(op);
		auto & so = op_cont->op;
		auto & grid = so.grid();

		for (auto i: range(nvals)) {
			len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
			len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
			len_t ck = static_cast<len_t>(coords[i].k - grid.is(2) + 2);
			// boxmg likes positive stencil coefficients
			if (coords[i].dir != 0) vals[i] = -1*vals[i];
			// if (coords[i].dir == 0)
			// log::error << ci << " " << cj << " -> " << coords[i].dir << " => " << vals[i]<< std::endl;
			so(static_cast<len_t>(ci),
			   static_cast<len_t>(cj),
			   static_cast<len_t>(ck),
			   static_cast<xxvii_pt>(coords[i].dir)) = vals[i];
		}
	}


	void bmg3_operator_apply(bmg3_operator op, const double *x, double *b)
	{
		using namespace cedar::cdr3;

		auto * op_cont = reinterpret_cast<op_container*>(op);
		auto & sop = op_cont->op;
		auto & xgf = op_cont->xgf;
		auto & bgf = op_cont->bgf;
		auto kman = op_cont->kman;

		int idx = 0;
		for (auto k : xgf.range(2)) {
			for (auto j : xgf.range(1)) {
				for (auto i : xgf.range(0)) {
					xgf(i,j,k) = x[idx];
					idx++;
				}
			}
		}

		kman->run<mpi::halo_exchange>(xgf);

		kman->run<mpi::matvec>(sop, xgf, bgf);

		idx = 0;
		for (auto k : bgf.range(2)) {
			for (auto j : bgf.range(1)) {
				for (auto i : bgf.range(0)) {
					b[idx] = bgf(i,j,k);
					idx++;
				}
			}
		}
	}


	void bmg3_operator_dump(bmg3_operator op)
	{
		using namespace cedar::cdr3;
		std::ofstream ofile;

		auto * op_cont = reinterpret_cast<op_container*>(op);
		auto & grid = op_cont->op.grid();

		ofile.open("op" + std::to_string(grid.coord(0)) + "-" + std::to_string(grid.coord(1)) + "-" + std::to_string(grid.coord(2)) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
		ofile << op_cont->op;
		ofile.close();
	}


	void bmg3_operator_destroy(bmg3_operator op)
	{
		using namespace cedar::cdr3;

		delete reinterpret_cast<op_container*>(op);
	}
}
