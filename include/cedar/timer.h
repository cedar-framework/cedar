#ifndef CEDAR_TIMER_H
#define CEDAR_TIMER_H

#include <cedar/global_manager.h>

namespace cedar {

inline void timer_begin(const std::string & label) { gman.get<gmant::timer>().begin(label); }
inline void timer_end(const std::string & label) { gman.get<gmant::timer>().end(label); }
inline void timer_up() { gman.get<gmant::timer>().up(); }
inline void timer_down() { gman.get<gmant::timer>().down(); }
inline void timer_pause() { gman.get<gmant::timer>().pause(); }
inline void timer_play() { gman.get<gmant::timer>().play(); }
inline void timer_redist(redist_comms comm) { gman.get<gmant::timer>().redist(comm); }
inline void timer_init(MPI_Comm comm) { gman.get<gmant::timer>().init(comm); }
inline void timer_save(const std::string & fname) { gman.get<gmant::timer>().save(fname); }


}

#endif
