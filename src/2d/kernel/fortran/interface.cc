#include "boxmg-common.h"

extern "C" {

	void print_error(char *string)
	{
		boxmg::log::error << string << std::endl;
	}

}
