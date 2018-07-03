#include <stdbool.h>

#include <cedar/types.h>
#include <cedar/2d/util/topo.h>
#include <cedar/config.h>

#include <cedar/2d/interface/c/operator.h>
#include <cedar/2d/interface/c/types.h>

extern "C"
{
	bmg2_operator bmg2_operator_create(bmg2_topo topo)
	{
		using namespace cedar;
		using namespace cedar::cdr2;

		auto grid = *(reinterpret_cast<std::shared_ptr<cedar::grid_topo>*>(topo));

		config conf("config.json");

		auto * op_cont = new op_container(grid, conf);

		return reinterpret_cast<bmg2_operator>(op_cont);
	}


	void bmg2_operator_set(bmg2_operator op, unsigned int nvals, grid_coord_2d coords[], double vals[])
	{
		using namespace cedar;
		using namespace cedar::cdr2;

		auto & so = reinterpret_cast<op_container*>(op)->op;
		auto & grid = so.grid();

		for (auto i: range(nvals)) {
			len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
			len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
			bmg2_dir dir = coords[i].dir;
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
			// if (coords[i].dir == 0)
			// log::error << ci << " " << cj << " -> " << coords[i].dir << " => " << vals[i]<< std::endl;
		}
	}


	void bmg2_operator_apply(bmg2_operator op, const double *x, double *b)
	{
		using namespace cedar::cdr2;

		auto *op_cont = reinterpret_cast<op_container*>(op);
		auto & sop = op_cont->op;

		auto kman = op_cont->kman;
		auto & halo_service = kman->services().get<mpi::halo_exchange>();

		auto & xgf = op_cont->xgf;
		auto & bgf = op_cont->bgf;

		int idx = 0;
		for (auto j : xgf.range(1)) {
			for (auto i : xgf.range(0)) {
				xgf(i,j) = x[idx];
				idx++;
			}
		}

		halo_service.run(xgf);
		kman->run<mpi::matvec>(sop, xgf, bgf);

		idx = 0;
		for (auto j : bgf.range(1)) {
			for (auto i : bgf.range(0)) {
				b[idx] = bgf(i,j);
				idx++;
			}
		}
	}


	void bmg2_operator_dump(bmg2_operator op)
	{
		using namespace cedar::cdr2;
		std::ofstream ofile;

		auto * op_cont = reinterpret_cast<op_container*>(op);
		auto & grid = op_cont->op.grid();

		ofile.open("op" + std::to_string(grid.coord(0)) + "-" + std::to_string(grid.coord(1)) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
		ofile << op_cont->op;
		ofile.close();

	}


	void bmg2_operator_destroy(bmg2_operator op)
	{
		using namespace cedar::cdr2;
		auto * op_cont = reinterpret_cast<op_container*>(op);
		delete op_cont;
	}
}
