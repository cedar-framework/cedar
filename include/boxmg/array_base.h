#ifndef BOXMG_ARRAY_BASE_H
#define BOXMG_ARRAY_BASE_H

namespace boxmg {

template <typename len_type>
class array_base
{
public:
	virtual len_type len(int i) const = 0;
};

}

#endif
