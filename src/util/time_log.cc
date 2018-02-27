#include <cedar/util/log.h>
#include <cedar/util/time_log.h>

namespace cedar {
	time_log<machine_mode::MPI> tlog;
	time_log<machine_mode::SERIAL> tlog_ser;
	bool serial_timers = true;
}
