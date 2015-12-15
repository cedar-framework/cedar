#ifndef BOXMG_2D_CORE_GRID_STENCIL_H
#define BOXMG_2D_CORE_GRID_STENCIL_H

#include "boxmg/2d/inter/types.h"
#include "boxmg/2d/types.h"
#include "boxmg/array.h"
#include "boxmg/2d/boundary_iterator.h"

namespace boxmg { namespace bmg2d {

	class grid_stencil : public array<len_t, real_t, 3>
	{
	public:
		grid_stencil() {};
		grid_stencil(len_t nx, len_t ny, unsigned int nghosts=1, bool intergrid=false, bool symmetric=true, bool five_pt=false);

		using array<len_t,real_t,3>::index;
		using array<len_t,real_t,3>::operator();
		len_t index(len_t i, len_t j, dir direction) const
		{
			return this->index(i,j,static_cast<len_t>(direction));
		}

		len_t index(len_t i, len_t j, inter::dir direction) const
		{
			return this->index(i,j,static_cast<len_t>(direction));
		}


		real_t & operator()(std::array<len_t,2> idx, dir dir)
		{
			return (*this)(idx[0], idx[1], dir);
		}

		real_t & operator()(len_t i, len_t j, dir direction)
		{
			symmetric_hack(i,j,direction);
			#ifdef DEBUG
			return vec.at(index(i,j,direction));
			#else
			return vec.data()[index(i,j,direction)];
			#endif
		}

		real_t & operator()(len_t i, len_t j, inter::dir direction)
		{
			#ifdef DEBUG
			return vec.at(index(i,j,direction));
			#else
			return vec.data()[index(i,j,direction)];
			#endif
		}


		real_t operator()(len_t i, len_t j, inter::dir direction) const
		{
			#ifdef DEBUG
			return vec.at(index(i,j,direction));
			#else
			return vec.data()[index(i,j,direction)];
			#endif
		}

		real_t operator()(len_t i, len_t j, dir direction) const
		{
			symmetric_hack(i,j,direction);

			#ifdef DEBUG
			return vec.at(index(i,j,direction));
			#else
			return vec.data()[index(i,j,direction)];
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

		const virtual range_t<len_t> & range(int i) const;
		const virtual range_t<len_t> & grange(int i) const;
		virtual len_t shape(int i) const;

		/* AlignedVector<real_t> & vector() { return data_;} */
		/* const AlignedVector<real_t> & vector() const { return data_;} */

	private:
		bool symmetric;
		bool five_pt_;
		int num_ghosts;
		std::array<range_t<len_t>, 2> range_;
		std::array<range_t<len_t>, 2> grange_;
	};

}}

#endif
