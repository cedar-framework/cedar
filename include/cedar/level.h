#ifndef CEDAR_LEVEL_H
#define CEDAR_LEVEL_H

namespace cedar
{

template<class sten, class solver_types>
struct level
{
	template<class stencil>
	using stencil_op = typename solver_types::template stencil_op<stencil>;
	using grid_func = typename solver_types::grid_func;
	using prolong_op = typename solver_types::prolong_op;
	using restrict_op = typename solver_types::restrict_op;
	using relax_stencil = typename solver_types::relax_stencil;

level(stencil_op<sten> & A) : A(A) {}
	template<class... Args>
	level(Args&&... args) : Adata(std::forward<Args>(args)...),
		A(Adata), P(std::forward<Args>(args)...), x(std::forward<Args>(args)...),
		res(std::forward<Args>(args)...), b(std::forward<Args>(args)...) {}
	stencil_op<sten> Adata;
	stencil_op<sten> & A;
	prolong_op  P;
	restrict_op R;
	grid_func x;
	grid_func res;
	grid_func b;
	std::array<relax_stencil, 2> SOR;

	std::function<void(const stencil_op<sten> & A, grid_func & x, const grid_func & b)> presmoother;
	std::function<void(const stencil_op<sten> & A, grid_func & x, const grid_func & b)> postsmoother;
};

}

#endif
