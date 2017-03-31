#include <stdbool.h>

#include <cedar/types.h>
#include <cedar/3d/util/topo.h>
#include <cedar/3d/mpi/stencil_op.h>
#include <cedar/3d/kernel/mpi/registry.h>

#include <cedar/3d/interface/c/operator.h>

extern "C"
{
	cdr3_operator cdr3_operator_create(cdr3_topo topo)
	{
		using namespace cedar::cdr3;

		auto & grid = *(reinterpret_cast<std::shared_ptr<cedar::grid_topo>*>(topo));

		mpi::stencil_op *sop = new mpi::stencil_op(grid);
		auto & sten = sop->stencil();

		return reinterpret_cast<cdr3_operator>(sop);
	}


	cdr3_operator cdr3_operator_create_fivept(cdr3_topo topo)
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		auto & grid = *(reinterpret_cast<std::shared_ptr<grid_topo>*>(topo));

		mpi::stencil_op *sop = new mpi::stencil_op(grid);
		auto & sten = sop->stencil();
		sten.five_pt() = true;

		return reinterpret_cast<cdr3_operator>(sop);
	}


	void cdr3_operator_set(cdr3_operator op, unsigned int nvals, grid_coord_3d coords[], double vals[])
	{
		using namespace cedar;
		using namespace cedar::cdr3;

		auto & sten = reinterpret_cast<mpi::stencil_op*>(op)->stencil();
		auto & grid = reinterpret_cast<mpi::stencil_op*>(op)->grid();

		for (auto i: range(nvals)) {
			len_t ci = static_cast<len_t>(coords[i].i - grid.is(0) + 2);
			len_t cj = static_cast<len_t>(coords[i].j - grid.is(1) + 2);
			len_t ck = static_cast<len_t>(coords[i].k - grid.is(2) + 2);
			// boxmg likes positive stencil coefficients
			if (coords[i].dir != 0) vals[i] = -1*vals[i];
			// if (coords[i].dir == 0)
			// log::error << ci << " " << cj << " -> " << coords[i].dir << " => " << vals[i]<< std::endl;
			sten(static_cast<len_t>(ci),
			     static_cast<len_t>(cj),
			     static_cast<len_t>(ck),
			     static_cast<dir>(coords[i].dir)) = vals[i];
		}
	}


	void cdr3_operator_apply(cdr3_operator op, const double *x, double *b)
	{
		using namespace cedar::cdr3;

		auto *sop = reinterpret_cast<mpi::stencil_op*>(op);
		auto grid = sop->grid_ptr();

		mpi::grid_func xgf(grid);
		int idx = 0;
		for (auto k : xgf.range(2)) {
			for (auto j : xgf.range(1)) {
				for (auto i : xgf.range(0)) {
					xgf(i,j,k) = x[idx];
					idx++;
				}
			}
		}

		std::shared_ptr<kernel::mpi::registry> kreg;
		kreg = std::dynamic_pointer_cast<kernel::mpi::registry>(sop->get_registry());
		xgf.halo_ctx = sop->halo_ctx;
		kreg->halo_exchange(xgf);

		mpi::grid_func bgf(grid);
		sop->apply(xgf, bgf);

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


	void cdr3_operator_dump(cdr3_operator op)
	{
		using namespace cedar::cdr3;
		std::ofstream ofile;

		auto &grid = reinterpret_cast<mpi::stencil_op*>(op)->grid();

		ofile.open("op" + std::to_string(grid.coord(0)) + "-" + std::to_string(grid.coord(1)) + "-" + std::to_string(grid.coord(2)) + ".txt", std::ios::out | std::ios::trunc | std::ios::binary);
		ofile << *(reinterpret_cast<mpi::stencil_op*>(op));
		ofile.close();

	}


	void cdr3_operator_destroy(cdr3_operator op)
	{
		using namespace cedar::cdr3;

		delete reinterpret_cast<mpi::stencil_op*>(op);
	}
}
