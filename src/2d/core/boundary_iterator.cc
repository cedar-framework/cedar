#include "boundary_iterator.h"

#include "grid_stencil.h"

using namespace boxmg::bmg2d::core;

BoundaryIterator::BoundaryIterator(const GridStencil &s, Dir dir):
	begin_(s,dir), end_(s,dir)
{
	if (dir == Dir::S) {
		begin_.set(s.range(0).front(), s.range(1).front());
		end_.set(s.range(0).back()+1, s.range(1).front());
	} else if (dir == Dir::W) {
		begin_.set(s.range(0).front(), s.range(1).front());
		end_.set(s.range(0).front(), s.range(1).back()+1);
	} else if (dir == Dir::N) {
		begin_.set(s.range(0).front(), s.range(1).back());
		end_.set(s.range(0).back()+1, s.range(1).back());
	} else if (dir == Dir::E) {
		begin_.set(s.range(0).back(), s.range(1).front());
		end_.set(s.range(0).back(), s.range(1).back()+1);
	}
}


const BoundaryIterator::iterator & BoundaryIterator::iterator::operator ++ ()
{
	if (dir == Dir::S or dir == Dir::N) {
		idx_[0]++;
	} else {
		idx_[1]++;
	}
	return *this;
}


BoundaryIterator::iterator BoundaryIterator::iterator::operator ++(int)
{
	BoundaryIterator::iterator copy(*this);

	if (dir == Dir::S or dir == Dir::N) {
		idx_[0]++;
	} else {
		idx_[1]++;
	}

	return copy;
}
