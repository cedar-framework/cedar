#ifndef BOXMG_2D_CORE_GRID_STENCIL_H
#define BOXMG_2D_CORE_GRID_STENCIL_H

#include "inter/types.h"
#include "core/types.h"
#include "core/stencil.h"
#include "core/array.h"
#include "core/boundary_iterator.h"

namespace boxmg { namespace bmg2d {

	class GridStencil : public Array<len_t, real_t>
	{
	public:
		GridStencil() {};
		GridStencil(len_t nx, len_t ny, unsigned int nghosts=1, bool intergrid=false, bool symmetric=true, bool five_pt=false);

		using Array<len_t,real_t>::index;
		len_t index(len_t i, len_t j, Dir dir) const
		{
			return static_cast<int>(dir)*ostride + i*stride_[0] + j*stride_[1];
		}

		len_t index(len_t i, len_t j, inter::Dir dir) const
		{
			return static_cast<int>(dir)*ostride + i*stride_[0] + j*stride_[1];
		}


		real_t & operator()(std::array<len_t,2> idx, Dir dir)
		{
			return (*this)(idx[0], idx[1], dir);
		}

		real_t & operator()(len_t i, len_t j, Dir dir)
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			#endif
			symmetric_hack(i,j,dir);
			#ifdef DEBUG
			return data_.at(index(i,j,dir));
			#else
			return data_.data()[index(i,j,dir)];
			#endif
		}

		real_t & operator()(len_t i, len_t j, inter::Dir dir)
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			return data_.at(index(i,j,dir));
			#else
			return data_.data()[index(i,j,dir)];
			#endif
		}


		real_t operator()(len_t i, len_t j, inter::Dir dir) const
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			return data_.at(index(i,j,dir));
			#else
			return data_.data()[index(i,j,dir)];
			#endif
		}

		real_t operator()(len_t i, len_t j, Dir dir) const
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			#endif

			symmetric_hack(i,j,dir);

			#ifdef DEBUG
			return data_.at(index(i,j,dir));
			#else
			return data_.data()[index(i,j,dir)];
			#endif
		}


		const virtual BoundaryIterator boarder(Dir dir) const
		{
			BoundaryIterator ret(*this, dir);
			return ret;
		}

		void symmetric_hack(len_t &i, len_t &j, Dir &dir) const
		{
			if (symmetric and (dir == Dir::N or dir == Dir::SE or
							   dir == Dir::NE or dir == Dir::E or
					           dir == Dir::NW)) {
				if (dir == Dir::SE) { i++; dir = Dir::NW; }
				else if (dir == Dir::N) { j++; dir = Dir::S; }
				else if (dir == Dir::NE) { i++; j++; dir = Dir::SW; }
				else if (dir == Dir::E) { i++; dir = Dir::W; }
				else if (dir == Dir::NW) {j++;}
			}
		}

		bool five_pt() const
		{
			return five_pt_;
		}

		bool & five_pt()
		{
			return five_pt_;
		}

		len_t & stride(int dim)
		{
			if (dim < 2) return stride_[dim];
			else return ostride;
		}

		len_t stride(int dim) const
		{
			if (dim < 2) return stride_[dim];
			else return ostride;
		}

		/* AlignedVector<real_t> & vector() { return data_;} */
		/* const AlignedVector<real_t> & vector() const { return data_;} */

	private:
		bool symmetric;
		bool five_pt_;
		len_t ostride;
	};

}}

#endif
