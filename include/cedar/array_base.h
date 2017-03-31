#ifndef CEDAR_ARRAY_BASE_H
#define CEDAR_ARRAY_BASE_H

namespace cedar {

template <typename len_type>
class array_base
{
public:
	virtual len_type len(unsigned short i) const = 0;
};

}

#endif
