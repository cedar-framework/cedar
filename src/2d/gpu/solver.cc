#include <cedar/2d/gpu/solver.h>

namespace cedar::cdr2::gpu::mpi {

    template<class fsten>
    std::size_t solver<fsten>::compute_num_levels(stencil_op<fsten> & fop)
    {
        int ng;

        kman->template run<setup_nog>(fop.grid(), this->settings.min_coarse, &ng);

        return ng;
    }
    template std::size_t solver<five_pt>::compute_num_levels(stencil_op<five_pt> & fop);
    template std::size_t solver<nine_pt>::compute_num_levels(stencil_op<nine_pt> & fop);


    template<class fsten>
    void solver<fsten>::solve(const cdr2::gpu::mpi::grid_func & b, cdr2::gpu::mpi::grid_func & x)
    {
        auto & bd = const_cast<cdr2::gpu::mpi::grid_func&>(b);
        auto & halo_service = kman->services().template get<halo_exchange>();
        halo_service.run(bd);
        return parent::solve(b, x);
    }
    template void solver<five_pt>::solve(const cdr2::gpu::mpi::grid_func & b, cdr2::gpu::mpi::grid_func & x);
    template void solver<nine_pt>::solve(const cdr2::gpu::mpi::grid_func & b, cdr2::gpu::mpi::grid_func & x);


    template<class fsten>
    grid_func solver<fsten>::solve(const grid_func & b)
    {
        auto & bd = const_cast<cdr2::gpu::mpi::grid_func&>(b);
        auto & halo_service = kman->services().template get<halo_exchange>();
        halo_service.run(bd);
        return parent::solve(b);
    }
    template grid_func solver<five_pt>::solve(const grid_func & b);
    template grid_func solver<nine_pt>::solve(const grid_func & b);

    /* Just use the default multilevel coarse solve */
    // template<class fsten>
    // void solver<fsten>::setup_cg_solve()
    // {
    //     auto params = build_kernel_params(*(this->conf));
    //     auto & cop = this->levels.get(this->levels.size() - 1).A;

    //     auto cg_conf = this->settings.coarse_config;

    //     // std::array<int, 2> choice{{1,1}};

    //     kman->template setup<solve_cg>(cop, ABD);
    //     auto kernels = get_kernels();
    //     coarse_solver = [&, kernels](grid_func &x, const grid_func & b) {
    //         kernels->template run<solve_cg>(x, b, ABD, bbd);
    //     };

    //     // log::status << "Redistributing to " << choice[0] << " x " << choice[1] << " cores" << std::endl;
    //     // if (this->settings.coarse_solver == ml_settings::cg_type::lu) {
    //     //     this->coarse_solver = create_redist_solver<cholesky_solver>(kman,
    //     //                                                                 *this->conf,
    //     //                                                                 cop,
    //     //                                                                 cg_conf,
    //     //                                                                 choice);
    //     // } else {
    //     //     using inner_solver = multilevel_wrapper<cdr2::solver<nine_pt>>;
    //     //     this->coarse_solver = create_redist_solver<inner_solver>(kman,
    //     //                                                              *this->conf,
    //     //                                                              cop,
    //     //                                                              cg_conf,
    //     //                                                              choice);
    //     // }
    // }
 // template void solver<five_pt>::setup_cg_solve();
    // template void solver<nine_pt>::setup_cg_solve();
    template <class fsten>
    void solver<fsten>::setup_cg_solve() {
        auto& level = this->levels.get(this->levels.size() - 1);
        // auto& cop = this->levels.get(this->levels.size() - 1).A;
        // auto& copr = this->levels.get(this->levels.size() - 1).SOR[0];

        kman->template setup<point_relax>(level.A, level.SOR[0]);

        auto kernels = this->get_kernels();
        this->coarse_solver = [&, kernels](grid_func &x, const grid_func & b) {
            for (int i = 0; i < 1; ++i) {
                //kernels->template run<point_relax>(cop, x, b, copr, cycle::Dir::DOWN);
                kernels->template run<point_relax>(level.A, x, b, level.SOR[0], cycle::Dir::DOWN);
            }
        };
    }
    template void solver<five_pt>::setup_cg_solve();
    template void solver<nine_pt>::setup_cg_solve();

    template<class fsten>
    void solver<fsten>::setup_space(std::size_t nlevels)
    {
        service_manager<stypes> & sman = kman->services();
        this->levels.init(sman.get<mempool>(),
                          nlevels);
        for (auto i : range<std::size_t>(nlevels-1)) {

            auto & fgrid = this->get_grid(i);
            if (i == 0)
                fgrid.grow(nlevels);

            int kc = nlevels - i - 2;

            auto cgrid = std::make_shared<grid_topo>(fgrid.get_igrd(), kc, nlevels);
            cgrid->comm = fgrid.comm;

            len_t NLxg = fgrid.nlocal(0) - 2;
            len_t NLyg = fgrid.nlocal(1) - 2;
            len_t NGxg = (fgrid.nglobal(0) - 1) / 2 + 2;
            len_t NGyg = (fgrid.nglobal(1) - 1) / 2 + 2;

            cgrid->nglobal(0) = NGxg;
            cgrid->nglobal(1) = NGyg;

            if ((fgrid.is(0) % 2) == 1) {
                cgrid->is(0) = (fgrid.is(0) + 1) / 2;
                NLxg = (NLxg + 1) / 2;
            } else {
                cgrid->is(0) = fgrid.is(0)/2 + 1;
                if (NLxg % 2 == 1) NLxg = (NLxg-1)/2;
                else NLxg = (NLxg+1)/2;
            }


            if (fgrid.is(1) % 2 == 1) {
                cgrid->is(1) = (fgrid.is(1)+1) / 2;
                NLyg = (NLyg+1) / 2;
            } else {
                cgrid->is(1) = fgrid.is(1) / 2 + 1;
                if (NLyg % 2 == 1) NLyg = (NLyg - 1) / 2;
                else NLyg = (NLyg+1)/2;
            }

            cgrid->nlocal(0) = NLxg + 2;
            cgrid->nlocal(1) = NLyg + 2;

            cgrid->nproc(0) = fgrid.nproc(0);
            cgrid->nproc(1) = fgrid.nproc(1);
            cgrid->coord(0) = fgrid.coord(0);
            cgrid->coord(1) = fgrid.coord(1);

            this->levels.add(sman.get<mempool>(),
                             cgrid);
        }
        setup_halo();
        {
            auto & sop = this->levels.template get<fsten>(0).A;
            auto & halo_service = kman->services().template get<halo_exchange>();
            halo_service.run(sop);
        }

        /* Allocate GPU memory for everything if we have a device */
        setup_gpu();
        for (std::size_t i = 0; i < this->nlevels(); i++) {
            if (i == 0) {
                auto& level = this->levels.template get<fsten>(i);
                level.A.ensure_gpu();
                level.P.ensure_gpu();
                level.x.ensure_gpu();
                level.res.ensure_gpu();
                level.b.ensure_gpu();
                auto resb = level.res.to_flat_buffer();
                std::cerr << level.res.data() << std::endl;
                std::cerr << "res: " << resb.tostr_unformatted() << std::endl;
            } else {
                auto& level = this->levels.get(i);
                level.A.ensure_gpu();
                level.P.ensure_gpu();
                level.x.ensure_gpu();
                level.res.ensure_gpu();
                level.b.ensure_gpu();
                auto resb = level.res.to_flat_buffer();
                std::cerr << level.res.data() << std::endl;
                std::cerr << "res: " << resb.tostr_unformatted() << std::endl;
            }
        }
    }
    template void solver<five_pt>::setup_space(std::size_t nlevels);
    template void solver<nine_pt>::setup_space(std::size_t nlevels);

    template<class fsten>
    void solver<fsten>::setup_halo()
    {
        auto & sop = this->levels.template get<fsten>(0).A;

        std::vector<topo_ptr> topos;
        topos.push_back(sop.grid_ptr());

        for (std::size_t i = 1; i < this->nlevels(); i++)
            topos.push_back(this->levels.get(i).A.grid_ptr());

        service_manager<stypes> & sman = kman->services();
        auto & halo_service = sman.get<halo_exchange>();
        halo_service.setup(topos);
    }
    template void solver<five_pt>::setup_halo();
    template void solver<nine_pt>::setup_halo();

    template <class fsten>
    void solver<fsten>::setup_gpu() {
        ftl::device::autodetect();
        log::status << "Running on " << ftl::device::get_name() << std::endl;

        /* Load kernels */
        if (!ftl::load_kernels(true)) {
            throw std::runtime_error("FTL: Failed to load kernels.");
        }
        log::status << "Loaded and compiled all GPU kernels" << std::endl;
    }
}
