#ifndef BOXMG_TIMER_H
#define BOXMG_TIMER_H

#include <string>

namespace boxmg {

	class timer
	{
	public:
		template<class T>
			timer(T&& name): name(std::forward<T>(name)) {};
		void begin();
		void end();
		double time();
	private:
		std::string name;
		double starttime, endtime;
	};

}

#endif
