#ifndef CEDAR_3D_RELAX_PLANES_H
#define CEDAR_3D_RELAX_PLANES_H

#include <cedar/types.h>
#include <cedar/3d/types.h>
#include <cedar/kernels/plane_relax.h>
#include <cedar/2d/solver.h>


namespace cedar { namespace cdr3 {


template<relax_dir rdir>
struct relax_dir_name {static constexpr const char* value = "";};

template<>
struct relax_dir_name<relax_dir::xy> {static constexpr const char* value = "xy";};

template<>
struct relax_dir_name<relax_dir::xz> {static constexpr const char* value = "xz";};

template<>
struct relax_dir_name<relax_dir::yz> {static constexpr const char* value = "yz";};

template<class fsten>
using slv2_ptr = std::unique_ptr<cdr2::solver<fsten>>;


int log_begin(bool log_planes, int ipl, const std::string & suff);
void log_end(bool log_planes, int ipl, int lvl);

template<relax_dir rdir, class sten3>
void copy_rhs(const stencil_op<sten3> & so,
              const grid_func & x,
              const grid_func & b,
              cdr2::grid_func & b2,
              int ipl);

template<relax_dir rdir>
void copy32(const grid_func & x, cdr2::grid_func & x2, int ipl);

template<relax_dir rdir>
void copy23(const cdr2::grid_func & x2, grid_func &x, int ipl);


template<class sten3>
using sten3_to_sten2 = typename std::conditional<std::is_same<sten3, seven_pt>::value,
                                                 cdr2::five_pt, cdr2::nine_pt>::type;
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
                cdr2::stencil_op<sten2> & so2);

template<relax_dir rdir>
void copy_coeff(const stencil_op<seven_pt> & so3,
                cdr2::stencil_op<cdr2::five_pt> & so2)
{
	using namespace cdr2;

	// TODO: reorder loops for performance
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
                cdr2::stencil_op<cdr2::nine_pt> & so2)
{
	using namespace cdr2;

	// TODO: reorder loops for performance
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
		if (rdir == relax_dir::xy) {
			for (auto k : so.range(2)) {
				auto so2_ptr = std::make_unique<cdr2::stencil_op<sten2>>(so.shape(0), so.shape(1));
				auto & so2 = *so2_ptr;
				copy_coeff<rdir>(so, so2);
				auto conf2 = this->params->plane_config;

				planes.emplace_back(std::make_unique<cdr2::solver<sten2>>(so2, conf2));
				planes.back()->give_op(std::move(so2_ptr));
				planes.back()->levels.template get<sten2>(0).x = cdr2::grid_func(so.shape(0), so.shape(1));
				planes.back()->levels.template get<sten2>(0).b = cdr2::grid_func(so.shape(0), so.shape(1));
			}
		} else if (rdir == relax_dir::xz) {
			for (auto j : so.range(1)) {
				auto so2_ptr = std::make_unique<cdr2::stencil_op<sten2>>(so.shape(0), so.shape(2));
				auto & so2 = *so2_ptr;
				copy_coeff<rdir>(so, so2);
				auto conf2 = this->params->plane_config;

				planes.emplace_back(std::make_unique<cdr2::solver<sten2>>(so2, conf2));
				planes.back()->give_op(std::move(so2_ptr));
				planes.back()->levels.template get<sten2>(0).x = cdr2::grid_func(so.shape(0), so.shape(2));
				planes.back()->levels.template get<sten2>(0).b = cdr2::grid_func(so.shape(0), so.shape(2));
			}
		} else if (rdir == relax_dir::yz) {
			for (auto i : so.range(0)) {
				auto so2_ptr = std::make_unique<cdr2::stencil_op<sten2>>(so.shape(1), so.shape(2));
				auto & so2 = *so2_ptr;
				copy_coeff<rdir>(so, so2);
				auto conf2 = this->params->plane_config;

				planes.emplace_back(std::make_unique<cdr2::solver<sten2>>(so2, conf2));
				planes.back()->give_op(std::move(so2_ptr));
				planes.back()->levels.template get<sten2>(0).x = cdr2::grid_func(so.shape(1), so.shape(2));
				planes.back()->levels.template get<sten2>(0).b = cdr2::grid_func(so.shape(1), so.shape(2));
			}
		} else {
			log::error << "invalid relax_dir for plane relaxation" << std::endl;
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
};

}}

#endif
