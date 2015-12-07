#include <boxmg/2d/boundary_iterator.h>
#include <boxmg/2d/grid_stencil.h>


using namespace boxmg::bmg2d;

BoundaryIterator::BoundaryIterator(const grid_stencil &s, dir direction):
	begin_(s,direction), end_(s,direction)
{
	if (direction == dir::S) {
		begin_.set(s.range(0).front(), s.range(1).front());
		end_.set(s.range(0).back()+1, s.range(1).front());
	} else if (direction == dir::W) {
		begin_.set(s.range(0).front(), s.range(1).front());
		end_.set(s.range(0).front(), s.range(1).back()+1);
	} else if (direction == dir::N) {
		begin_.set(s.range(0).front(), s.range(1).back());
		end_.set(s.range(0).back()+1, s.range(1).back());
	} else if (direction == dir::E) {
		begin_.set(s.range(0).back(), s.range(1).front());
		end_.set(s.range(0).back(), s.range(1).back()+1);
	}
}


const BoundaryIterator::iterator & BoundaryIterator::iterator::operator ++ ()
{
	if (direction == dir::S or direction == dir::N) {
		idx_[0]++;
	} else {
		idx_[1]++;
	}
	return *this;
}


BoundaryIterator::iterator BoundaryIterator::iterator::operator ++(int)
{
	BoundaryIterator::iterator copy(*this);

	if (direction == dir::S or direction == dir::N) {
		idx_[0]++;
	} else {
		idx_[1]++;
	}

	return copy;
}
