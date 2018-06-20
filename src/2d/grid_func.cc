#include <algorithm>
#include <random>

#include <cedar/2d/grid_func.h>

using namespace cedar;
using namespace cedar::cdr2;


grid_func & grid_func::operator=(grid_func &&gf)
{
	vec = std::move(gf.vec);
	strides = std::move(gf.strides);
	extents = std::move(gf.extents);
	num_ghosts = gf.num_ghosts;
	range_ = std::move(gf.range_);
	grange_ = std::move(gf.grange_);

	return *this;
}


grid_func::grid_func(len_t nx, len_t ny, unsigned int nghosts) :
	array<real_t,2>(nx+2*nghosts, ny+2*nghosts)
{
	num_ghosts = nghosts;
	range_[0] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));

	grange_[0] = cedar::range(static_cast<len_t>(0), nx + 2*nghosts);
	grange_[1] = cedar::range(static_cast<len_t>(0), ny + 2*nghosts);
}


grid_func::grid_func(real_t *ext_data, len_t nx, len_t ny, unsigned int nghosts) :
	array<real_t,2>(ext_data, nx+2*nghosts, ny+2*nghosts)
{
	num_ghosts = nghosts;
	range_[0] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(nx + nghosts));
	range_[1] = cedar::range(static_cast<len_t>(nghosts), static_cast<len_t>(ny + nghosts));

	grange_[0] = cedar::range(static_cast<len_t>(0), nx + 2*nghosts);
	grange_[1] = cedar::range(static_cast<len_t>(0), ny + 2*nghosts);
}


grid_func grid_func::ones(len_t nx, len_t ny)
{
	grid_func ret(nx,ny);

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


grid_func grid_func::zeros(len_t nx, len_t ny)
{
	grid_func ret(nx,ny);

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


grid_func grid_func::random(len_t nx, len_t ny)
{
	grid_func ret(nx,ny);

	std::uniform_real_distribution<double> unif(0, 1);
	std::random_device rand_dev;

	std::default_random_engine rand_engine(rand_dev());

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = unif(rand_engine);
		}
	}

	return ret;
}


grid_func grid_func::like(const grid_func &likeable)
{
	grid_func ret(likeable.shape(0), likeable.shape(1));

	return ret;
}


grid_func grid_func::zeros_like(const grid_func &like)
{
	grid_func ret(like.shape(0), like.shape(1));

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 0.0;
		}
	}

	return ret;
}


grid_func grid_func::ones_like(const grid_func &like)
{
	grid_func ret(like.shape(0), like.shape(1));

	for (auto j: ret.grange(1)) {
		for (auto i: ret.grange(0)) {
			ret(i,j) = 1.0;
		}
	}

	return ret;
}


cedar::real_t grid_func::inf_norm() const
{
	auto abs_compare = [](real_t a, real_t b){
		return (std::abs(a) < std::abs(b));
	};

	real_t cmax = 0;

	for (auto j : this->range(1)) {
		for (auto i : this->range(0)) {
			if (abs_compare(cmax, (*this)(i,j)))
				cmax = (*this)(i,j);
		}
	}

	return cmax;
}


grid_func & grid_func::operator-=(const grid_func &rhs)
{
	auto jj = rhs.range(1).begin();
	for (auto j: this->range(1)) {
		auto ii = rhs.range(0).begin();
		for (auto i: this->range(0)) {
			(*this)(i,j) -= rhs(*ii,*jj);
			++ii;
		}
		++jj;
	}

	return *this;
}


namespace cedar { namespace cdr2 {
std::ostream & operator<<(std::ostream &os, const grid_func & obj)
{
	unsigned int width = 4;

	os << std::setprecision(7);


	for (auto j: obj.range(1)) {
		for (auto i: obj.range(0)) {
			os << std::setw(width) << i+1 << ", " << std::setw(width) << j+1 << ", "
			   << std::scientific << obj(i,j) << '\n';
		}
	}
	return os;
}
}}
