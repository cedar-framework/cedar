#ifndef BOXMG_2D_CORE_GRID_STENCIL_H
#define BOXMG_2D_CORE_GRID_STENCIL_H

#include "inter/types.h"
#include "core/types.h"
#include "core/stencil.h"
#include "core/array.h"
#include "core/boundary_iterator.h"

namespace boxmg { namespace bmg2d {

	class grid_stencil : public Array<len_t, real_t>
	{
	public:
		grid_stencil() {};
		grid_stencil(len_t nx, len_t ny, unsigned int nghosts=1, bool intergrid=false, bool symmetric=true, bool five_pt=false);

		using Array<len_t,real_t>::index;
		len_t index(len_t i, len_t j, dir direction) const
		{
			return static_cast<int>(direction)*ostride + i*stride_[0] + j*stride_[1];
		}

		len_t index(len_t i, len_t j, inter::dir direction) const
		{
			return static_cast<int>(direction)*ostride + i*stride_[0] + j*stride_[1];
		}


		real_t & operator()(std::array<len_t,2> idx, dir dir)
		{
			return (*this)(idx[0], idx[1], dir);
		}

		real_t & operator()(len_t i, len_t j, dir direction)
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			#endif
			symmetric_hack(i,j,direction);
			#ifdef DEBUG
			return data_.at(index(i,j,direction));
			#else
			return data_.data()[index(i,j,direction)];
			#endif
		}

		real_t & operator()(len_t i, len_t j, inter::dir direction)
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			return data_.at(index(i,j,direction));
			#else
			return data_.data()[index(i,j,direction)];
			#endif
		}


		real_t operator()(len_t i, len_t j, inter::dir direction) const
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			return data_.at(index(i,j,direction));
			#else
			return data_.data()[index(i,j,direction)];
			#endif
		}

		real_t operator()(len_t i, len_t j, dir direction) const
		{
			#ifdef DEBUG
			this->check_bounds(i,j);
			#endif

			symmetric_hack(i,j,direction);

			#ifdef DEBUG
			return data_.at(index(i,j,direction));
			#else
			return data_.data()[index(i,j,direction)];
			#endif
		}


		const virtual BoundaryIterator boarder(dir direction) const
		{
			BoundaryIterator ret(*this, direction);
			return ret;
		}

		void symmetric_hack(len_t &i, len_t &j, dir &direction) const
		{
			if (symmetric and (direction == dir::N or direction == dir::SE or
							   direction == dir::NE or direction == dir::E or
					           direction == dir::NW)) {
				if (direction == dir::SE) { i++; direction = dir::NW; }
				else if (direction == dir::N) { j++; direction = dir::S; }
				else if (direction == dir::NE) { i++; j++; direction = dir::SW; }
				else if (direction == dir::E) { i++; direction = dir::W; }
				else if (direction == dir::NW) {j++;}
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
