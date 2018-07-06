#ifndef CEDAR_3D_MPI_RELAX_PLANES_H
#define CEDAR_3D_MPI_RELAX_PLANES_H

#include <cedar/types.h>
#include <cedar/kernels/plane_relax.h>
#include <cedar/3d/mpi/types.h>
#include <cedar/2d/mpi/solver.h>
#include <cedar/3d/mpi/kernel_manager.h>

namespace cedar { namespace cdr3 { namespace mpi {

int log_begin(bool log_planes, int ipl, const std::string & suff);
void log_end(bool log_planes, int ipl, int lvl);

cdr2::mpi::kman_ptr master_kman(config & conf, int nplanes);
cdr2::mpi::kman_ptr worker_kman(int worker_id, config & conf, cdr2::mpi::kman_ptr master);

template<class fsten>
using slv2_ptr = std::unique_ptr<cdr2::mpi::solver<fsten>>;

template<relax_dir rdir, class sten3>
void copy_rhs(const stencil_op<sten3> & so,
              const grid_func & x,
              const grid_func & b,
              cdr2::mpi::grid_func & b2,
              int ipl);

template<relax_dir rdir>
void copy32(const grid_func & x, cdr2::mpi::grid_func & x2, int ipl);

template<relax_dir rdir>
void copy23(const cdr2::mpi::grid_func & x2, grid_func &x, int ipl);


template<relax_dir rdir, class sten3, class sten2>
void relax_planes(const stencil_op<sten3> & so, grid_func & x,
                  const grid_func & b, cycle::Dir cdir,
                  std::vector<slv2_ptr<sten2>> & planes)
{
	auto & conf = planes[0]->get_config();
	auto log_planes = conf.template get<bool>("log-planes", true);
	int lstart, lend, lstride;

	if (cdir == cycle::Dir::DOWN) {
		lstart = 1;
		lend = 3;
		lstride = 1;
	} else {
		lstart = 2;
		lend = 0;
		lstride = -1;
	}

	// red-black relaxation
	for (auto ipl_beg : range<int>(lstart, lend, lstride)) {
		for (auto ipl : range<int>(ipl_beg, planes.size() + 1, 2)) {
			auto & x2 = planes[ipl-1]->levels.template get<sten2>(0).x;
			auto & b2 = planes[ipl-1]->levels.template get<sten2>(0).b;

			copy32<rdir>(x, x2, ipl);
			copy_rhs<rdir>(so, x, b, b2, ipl);

			auto tmp = log_begin(log_planes, ipl, relax_dir_name<rdir>::value);
			planes[ipl-1]->solve(b2, x2);
			log_end(log_planes, ipl, tmp);

			copy23<rdir>(x2, x, ipl);
		}
	}
}


template<relax_dir rdir, class sten3, class sten2>
void copy_coeff(const stencil_op<sten3> & so3,
                cdr2::mpi::stencil_op<sten2> & so2);

template<relax_dir rdir>
void copy_coeff(const stencil_op<seven_pt> & so3,
                cdr2::mpi::stencil_op<cdr2::five_pt> & so2)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto k : so3.range(2)) {
			for (auto j : so3.grange(1)) {
				for (auto i : so3.grange(0)) {
					so2(i,j,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(i,j,five_pt::w) = so3(i,j,k,seven_pt::pw);
					so2(i,j,five_pt::s) = so3(i,j,k,seven_pt::ps);
				}
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto j : so3.range(1)) {
			for (auto k : so3.grange(2)) {
				for (auto i : so3.grange(0)) {
					so2(i,k,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(i,k,five_pt::w) = so3(i,j,k,seven_pt::pw);
					so2(i,k,five_pt::s) = so3(i,j,k,seven_pt::b);
				}
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto i : so3.range(0)) {
			for (auto k : so3.grange(2)) {
				for (auto j : so3.grange(1)) {
					so2(j,k,five_pt::c) = so3(i,j,k,seven_pt::p);
					so2(j,k,five_pt::w) = so3(i,j,k,seven_pt::ps);
					so2(j,k,five_pt::s) = so3(i,j,k,seven_pt::b);
				}
			}
		}
	}
}

template<relax_dir rdir>
void copy_coeff(const stencil_op<xxvii_pt> & so3,
                cdr2::mpi::stencil_op<cdr2::nine_pt> & so2)
{
	using namespace cdr2;

	if (rdir == relax_dir::xy) {
		for (auto k : so3.range(2)) {
			for (auto j : so3.grange(1)) {
				for (auto i : so3.grange(0)) {
					so2(i,j,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(i,j,nine_pt::w) = so3(i,j,k,xxvii_pt::pw);
					so2(i,j,nine_pt::s) = so3(i,j,k,xxvii_pt::ps);
					so2(i,j,nine_pt::sw) = so3(i,j,k,xxvii_pt::psw);
					so2(i,j,nine_pt::nw) = so3(i,j,k,xxvii_pt::pnw);
				}
			}
		}
	} else if (rdir == relax_dir::xz) {
		for (auto j : so3.range(1)) {
			for (auto k : so3.grange(2)) {
				for (auto i : so3.grange(0)) {
					so2(i,k,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(i,k,nine_pt::w) = so3(i,j,k,xxvii_pt::pw);
					so2(i,k,nine_pt::s) = so3(i,j,k,xxvii_pt::b);
					so2(i,k,nine_pt::sw) = so3(i,j,k,xxvii_pt::bw);
					so2(i,k,nine_pt::nw) = so3(i,j,k,xxvii_pt::be);
				}
			}
		}
	} else if (rdir == relax_dir::yz) {
		for (auto i : so3.range(0)) {
			for (auto k : so3.grange(2)) {
				for (auto j : so3.grange(1)) {
					so2(j,k,nine_pt::c) = so3(i,j,k,xxvii_pt::p);
					so2(j,k,nine_pt::w) = so3(i,j,k,xxvii_pt::ps);
					so2(j,k,nine_pt::s) = so3(i,j,k,xxvii_pt::b);
					so2(j,k,nine_pt::sw) = so3(i,j,k,xxvii_pt::bs);
					so2(j,k,nine_pt::nw) = so3(i,j,k,xxvii_pt::bn);
				}
			}
		}
	}
}

template<relax_dir rdir>
class planes : public kernels::plane_relax<stypes, rdir>
{
public:
	void setup(const stencil_op<seven_pt> & so) override
	{
		this->setup_impl(so, fine_planes);
	}
	void setup(const stencil_op<xxvii_pt> & so) override
	{
		level_planes.emplace_back();
		level_map[so.shape(0)] = level_planes.size() - 1;
		this->setup_impl(so, level_planes.back());
	}
	void run(const stencil_op<seven_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }
	void run(const stencil_op<xxvii_pt> & so, grid_func & x,
	         const grid_func & b, cycle::Dir cycle_dir) override { this->run_impl(so, x, b, cycle_dir); }

	template<class sten3, class sten2>
	void setup_impl(const stencil_op<sten3> & so, std::vector<slv2_ptr<sten2>> & planes)
	{
		int nplanes = so.shape(2);
		auto rng = so.range(2);
		if (rdir == relax_dir::xz) {
			rng = so.range(1);
			nplanes = so.shape(1);
		} else if (rdir == relax_dir::yz) {
			rng = so.range(0);
			nplanes = so.shape(0);
		}
		auto topo2 = slice_topo(so.grid());
		auto conf2 = this->params->plane_config;
		auto log_planes = conf2->template get<bool>("log-planes", true);
		cdr2::mpi::kman_ptr master_kmans[2];
		master_kmans[0] = master_kman(*conf2, (nplanes / 2) + (nplanes % 2));
		master_kmans[1] = master_kman(*conf2, nplanes / 2);
		for (auto i : rng) {
			cdr2::mpi::kman_ptr kman2;
			if (i < 2)
				kman2 = master_kmans[i];
			else
				kman2 = worker_kman(i, *conf2, master_kmans[i % 2]);
			auto so2_ptr = std::make_unique<cdr2::mpi::stencil_op<sten2>>(topo2);
			auto & so2 = *so2_ptr;
			copy_coeff<rdir>(so, so2);

			auto tmp = log_begin(log_planes, i, relax_dir_name<rdir>::value);
			planes.emplace_back(std::make_unique<cdr2::mpi::solver<sten2>>(so2, conf2, kman2));
			log_end(log_planes, i, tmp);
			planes.back()->give_op(std::move(so2_ptr));
			planes.back()->levels.template get<sten2>(0).x = cdr2::mpi::grid_func(topo2);
			planes.back()->levels.template get<sten2>(0).b = cdr2::mpi::grid_func(topo2);
		}
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so, grid_func & x,
	              const grid_func & b, cycle::Dir cycle_dir)
	{

		if (std::is_same<sten, seven_pt>::value) {
			relax_planes<rdir>(so, x, b, cycle_dir, fine_planes);
		} else {
			len_t lsize = so.shape(0);
			relax_planes<rdir>(so, x, b, cycle_dir, level_planes[level_map[lsize]]);
		}
	}

protected:
	std::vector<std::vector<slv2_ptr<cdr2::nine_pt>>> level_planes;
	std::vector<slv2_ptr<cdr2::five_pt>> fine_planes;
	std::map<len_t, std::size_t> level_map;

	std::shared_ptr<grid_topo> slice_topo(const grid_topo & topo3)
	{
		auto igrd = std::make_shared<std::vector<len_t>>(NBMG_pIGRD);
		auto topo2 = std::make_shared<grid_topo>(igrd, 0, 1);

		topo2->nproc(2) = 1;
		if (rdir == relax_dir::xy) {
			MPI_Comm_split(topo3.comm, topo3.coord(2),
			               topo3.coord(1) * topo3.nproc(0) + topo3.coord(0), &topo2->comm);
			for (auto i : range<std::size_t>(2)) {
				topo2->nproc(i) = topo3.nproc(i);
				topo2->coord(i) = topo3.coord(i);
				topo2->is(i) = topo3.is(i);
				topo2->nlocal(i) = topo3.nlocal(i);
				topo2->nglobal(i) = topo3.nglobal(i);
			}

			auto & halo_service = this->services->template get<halo_exchange>();
			auto & dimx = halo_service.leveldims(0);
			auto & dimy = halo_service.leveldims(1);
			topo2->dimxfine.resize(topo2->nproc(0));
			topo2->dimyfine.resize(topo2->nproc(1));
			for (auto i : range<len_t>(topo2->nproc(0))) {
				topo2->dimxfine[i] = dimx(i, topo3.level());
			}

			for (auto j : range<len_t>(topo2->nproc(1))) {
				topo2->dimyfine[j] = dimy(j, topo3.level());
			}
		} else if (rdir == relax_dir::xz) {
			MPI_Comm_split(topo3.comm, topo3.coord(1),
			               topo3.coord(2) * topo3.nproc(0) + topo3.coord(0), &topo2->comm);
			for (auto i : range<std::size_t>(2)) {
				auto i3 = (i == 0) ? 0 : 2;
				topo2->nproc(i) = topo3.nproc(i3);
				topo2->coord(i) = topo3.coord(i3);
				topo2->is(i) = topo3.is(i3);
				topo2->nlocal(i) = topo3.nlocal(i3);
				topo2->nglobal(i) = topo3.nglobal(i3);
			}

			auto & halo_service = this->services->template get<halo_exchange>();
			auto & dimx = halo_service.leveldims(0);
			auto & dimy = halo_service.leveldims(2);
			topo2->dimxfine.resize(topo2->nproc(0));
			topo2->dimyfine.resize(topo2->nproc(1));
			for (auto i : range<len_t>(topo2->nproc(0))) {
				topo2->dimxfine[i] = dimx(i, topo3.level());
			}

			for (auto j : range<len_t>(topo2->nproc(1))) {
				topo2->dimyfine[j] = dimy(j, topo3.level());
			}
		} else if (rdir == relax_dir::yz) {
			MPI_Comm_split(topo3.comm, topo3.coord(0),
			               topo3.coord(2) * topo3.nproc(1) + topo3.coord(1), &topo2->comm);
			for (auto i : range<std::size_t>(2)) {
				auto i3 = i + 1;
				topo2->nproc(i) = topo3.nproc(i3);
				topo2->coord(i) = topo3.coord(i3);
				topo2->is(i) = topo3.is(i3);
				topo2->nlocal(i) = topo3.nlocal(i3);
				topo2->nglobal(i) = topo3.nglobal(i3);
			}

			auto & halo_service = this->services->template get<halo_exchange>();
			auto & dimx = halo_service.leveldims(1);
			auto & dimy = halo_service.leveldims(2);
			topo2->dimxfine.resize(topo2->nproc(0));
			topo2->dimyfine.resize(topo2->nproc(1));
			for (auto i : range<len_t>(topo2->nproc(0))) {
				topo2->dimxfine[i] = dimx(i, topo3.level());
			}

			for (auto j : range<len_t>(topo2->nproc(1))) {
				topo2->dimyfine[j] = dimy(j, topo3.level());
			}
		} else {
			log::error << "invalid relax_dir for planes" << std::endl;
		}

		return topo2;
	}
};

}}}

#endif
