#include <random>
#include <type_traits>
#include <limits>

#include <cuda_runtime.h>
#include <omp.h>

#include <cedar/kernels/residual.h>
#include <cedar/2d/types.h>
#include <cedar/2d/kernel_manager.h>
#include <cedar/util/timer_util.h>

#include "rescuda.h"

using namespace cedar;
using namespace cedar::cdr2;

extern "C" {
	using namespace cedar;
	void BMG2_SymStd_residual_intent(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
	                                 int*,int*,int*,int*,int*,int*,int*);
	void BMG2_SymStd_residual_offload(int*,real_t*,real_t*,real_t*,real_t*,len_t*,len_t*,
	                                  int*,int*,int*,int*,int*,int*,int*);

	void BMG_get_bc(int, int*);
}


class residual_intent : public kernels::residual<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
		this->run_impl(so, x, b, r);
	}


	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
		this->run_impl(so, x, b, r);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		int k = 0;
		int kf = 0;
		int ifd;
		int ibc;
		int nstncl = stencil_ndirs<sten>::value;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		auto & Ad = const_cast<stencil_op<sten>&>(so);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);

		BMG_get_bc(params->per_mask(), &ibc);
		BMG2_SymStd_residual_intent(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
		                            &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);
	}
};


class residual_omp : public kernels::residual<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
		this->run_impl(so, x, b, r);
	}


	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
		this->run_impl(so, x, b, r);
	}

	template<class sten>
	void run_impl(const stencil_op<sten> & so,
	              const grid_func & x,
	              const grid_func & b,
	              grid_func & r)
	{
		int k = 0;
		int kf = 0;
		int ifd;
		int ibc;
		int nstncl = stencil_ndirs<sten>::value;
		if (std::is_same<sten, five_pt>::value)
			ifd = 1;
		else
			ifd = 0;
		int irelax = 0;
		int irelax_sym = 0;
		int updown = 0;
		len_t ii = r.len(0);
		len_t jj = r.len(1);

		auto & Ad = const_cast<stencil_op<sten>&>(so);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);

		BMG_get_bc(params->per_mask(), &ibc);
		BMG2_SymStd_residual_offload(&k, Ad.data(), bd.data(), xd.data(), r.data(), &ii, &jj,
		                             &kf, &ifd, &nstncl, &ibc, &irelax, &irelax_sym, &updown);
	}
};


class residual_cuda : public kernels::residual<stypes>
{
	void run(const stencil_op<five_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
	}


	void run(const stencil_op<nine_pt> & so,
	         const grid_func & x,
	         const grid_func & b,
	         grid_func & r) override {
		int ilen = r.len(0);
		int jlen = r.len(1);
		auto & Ad = const_cast<stencil_op<nine_pt>&>(so);
		grid_func &xd = const_cast<grid_func&>(x);
		grid_func &bd = const_cast<grid_func&>(b);
		residual_cuda9(ilen, jlen, r.data(), Ad.data(), xd.data(), bd.data());
		cudaDeviceSynchronize();
	}
};


static void set_random(grid_func & x)
{
	using namespace cedar;

	std::mt19937 gen;
	gen.seed(0);
	std::uniform_real_distribution<real_t> dis;

	#pragma omp parallel for simd collapse(2)
	for (int j = 1; j < x.len(1); j++) {
		for (int i = 1; i < x.len(0); i++) {
			x(i,j) = dis(gen);
		}
	}
}


template<class sten>
static void set_random(stencil_op<sten> & so)
{
	using namespace cedar;

	std::mt19937 gen;
	gen.seed(0);
	std::uniform_real_distribution<real_t> dis;

	#pragma omp parallel for simd collapse(2)
	for (int j = 1; j < so.len(1); j++) {
		for (int i = 1; i < so.len(0); i++) {
			for (auto k : range(stencil_ndirs<sten>::value)) {
				so(i,j,static_cast<sten>(k)) = dis(gen);
			}
		}
	}
}


static void set_constant(grid_func & x, real_t val)
{
	using namespace cedar;

	len_t jlen = x.len(1);
	len_t ilen = x.len(0);
	real_t *data = x.data();
	#pragma omp target teams distribute parallel for simd collapse(2)
	for (int j = 1; j < jlen; j++) {
		for (int i = 1; i < ilen; i++) {
			data[j*ilen + i] = val;
		}
	}
}


template<class sten>
static void set_constant(stencil_op<sten> & so, real_t val)
{
	len_t jlen = so.len(1);
	len_t ilen = so.len(0);
	int klen = stencil_ndirs<sten>::value;
	real_t *data = so.data();
	#pragma omp target teams distribute parallel for simd collapse(2)
	for (int j = 1; j < jlen; j++) {
		for (int i = 1; i < ilen; i++) {
			for (int k = 0; k < klen; k++) {
				data[k*ilen*jlen + j*ilen + i] = val;
			}
		}
	}
}


template<class T>
struct managed_allocator
{
	static T * allocate(const std::size_t n)
	{
		void *pv;
		cudaMallocManaged(&pv, n*sizeof(T), cudaMemAttachGlobal);

		return static_cast<T*>(pv);
	}

	static void deallocate(T * const p, const std::size_t n)
	{
		cudaFree(p);
	}
};


template<class T>
struct explicit_allocator
{
	static T * allocate(const std::size_t n)
	{
		void *pv;
		cudaMalloc(&pv, n*sizeof(T));

		return static_cast<T*>(pv);
	}

	static void deallocate(T * const p, const std::size_t n)
	{
		cudaFree(p);
	}
};


template<class T>
struct default_allocator
{
	static T * allocate(const std::size_t n)
	{
		void *pv;

		posix_memalign(&pv, 64, n*sizeof(T));

		return static_cast<T*>(pv);
	}

	static void deallocate(T * const p, const std::size_t n)
	{
		free(p);
	}
};


template<class allocator>
grid_func allocate_vec(len_t nx, len_t ny)
{
	auto *addr = allocator::allocate((nx+2)*(ny+2));
	grid_func ret(addr, nx, ny);

	return ret;
}


template<class sten, class allocater>
stencil_op<sten> allocate_op(len_t nx, len_t ny)
{
	auto *addr = allocater::allocate((nx+2)*(ny+2)*stencil_ndirs<sten>::value);
	stencil_op<sten> ret(addr, nx, ny);

	return ret;
}


template<class sten>
std::function<stencil_op<sten>(len_t nx, len_t ny)> get_alloc_op(const std::string & memt)
{
	if (memt == "managed")
		return allocate_op<sten, managed_allocator<real_t>>;
	else if (memt == "explicit")
		return allocate_op<sten, explicit_allocator<real_t>>;
	else
		return allocate_op<sten, default_allocator<real_t>>;
}


std::function<grid_func(len_t nx, len_t ny)> get_alloc_vec(const std::string & memt)
{
	if (memt == "managed")
		return allocate_vec<managed_allocator<real_t>>;
	else if (memt == "explicit")
		return allocate_vec<explicit_allocator<real_t>>;
	else
		return allocate_vec<default_allocator<real_t>>;
}


template<class T>
void prefetch(T arr)
{
	int defdev = omp_get_default_device();
	cudaMemPrefetchAsync(arr.data(), arr.size()*sizeof(real_t), defdev, cudaStreamLegacy);
}


int main(int argc, char *argv[])
{
	config conf("restest.json");
	log::init(conf);
	auto ndofs = conf.getvec<len_t>("grid.n");
	auto nx = ndofs[0];
	auto ny = ndofs[1];

	auto kreg = build_kernel_manager(conf);
	kreg->add<residual, residual_intent>("intent");
	kreg->add<residual, residual_omp>("omp");
	kreg->add<residual, residual_cuda>("cuda");

	bool managed = conf.get<bool>("managed");
	auto kname = conf.get<std::string>("impl");
	std::string memt;
	if (managed)
		memt = "managed";
	else
		memt = "ats";
	if (kname == "cuda")
		memt = "explicit";
	auto alloc_op = get_alloc_op<nine_pt>(memt);
	auto alloc_vec = get_alloc_vec(memt);

	auto so = alloc_op(nx, ny);
	auto x = alloc_vec(nx, ny);
	auto b = alloc_vec(nx, ny);
	auto r = alloc_vec(nx, ny);

	kreg->set<residual>(kname);

	if (kname == "omp") {
		prefetch(so);
		prefetch(x);
		prefetch(b);
		prefetch(r);
		cudaDeviceSynchronize();
		set_constant(so, 4);
		set_constant(x, 2);
		set_constant(b, 1);
	} else if (kname == "cuda") {
		set_constant_cuda(so.data(), so.size(), 4);
		set_constant_cuda(x.data(), x.size(), 2);
		set_constant_cuda(b.data(), b.size(), 1);
		cudaDeviceSynchronize();
	} else {
		set_random(so);
		set_random(x);
		set_random(b);
	}

	int nruns = conf.get<int>("nruns");
	auto start = timer_util<machine_mode::SERIAL>::now();
	for (int i = 0; i < nruns; i++) {
		kreg->run<residual>(so, x, b, r);
	}
	auto end = timer_util<machine_mode::SERIAL>::now();
	auto elapsed = timer_util<machine_mode::SERIAL>::duration(start, end);
	elapsed = elapsed / nruns;
	std::cout << elapsed << " s " << nx * ny * 8 * (5+3) / 1e6 / elapsed << " MB/s" << std::endl;

	return 0;
}
