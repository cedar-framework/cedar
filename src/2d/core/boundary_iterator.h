#ifndef BOXMG_2D_CORE_BOUNDARY_ITERATOR
#define BOXMG_2D_CORE_BOUNDARY_ITERATOR

#include <array>

#include "boxmg-common.h"
#include "types.h"

namespace boxmg { namespace bmg2d {

class grid_stencil;

using i2d = std::array<len_t, 2>;
class BoundaryIterator {
public:
	class iterator {
		friend class BoundaryIterator;
	public:
		i2d operator *() const { return idx_; }
		const iterator &operator ++();
		iterator operator ++(int);

		bool operator ==(const iterator &other) const { return idx_[0] == other.idx_[0] and idx_[1] == other.idx_[1]; }
		bool operator !=(const iterator &other) const { return idx_[0] != other.idx_[0] or  idx_[1] != other.idx_[1]; }

		void set(len_t i, len_t j) { idx_[0] = i; idx_[1] = j; }

	protected:
	iterator(const grid_stencil & sten, Dir dir, i2d start): idx_(start), sten(sten), dir(dir) {}
	iterator(const grid_stencil &sten, Dir dir): sten(sten), dir(dir) {}

	private:
		i2d idx_;
		const grid_stencil & sten;
		Dir dir;
	};

	iterator begin() const { return begin_; }
	iterator end() const { return end_; }

	BoundaryIterator(const grid_stencil & s, Dir dir);

private:
	iterator begin_;
	iterator end_;
};

}}
#endif
