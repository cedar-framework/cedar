#include "types.h"
#include "log.h"

namespace boxmg { namespace log {

using lmap_t = std::map<std::string, unsigned int>;
std::unique_ptr<lmap_t> log_level = nullptr;
unsigned int level = 0;

LevelLogger memory("memory", Color::Modifier(Color::FG_BLUE));
LevelLogger status("status", Color::Modifier(Color::FG_DEFAULT));
LevelLogger info("info", Color::Modifier(Color::FG_DEFAULT));
LevelLogger error("error", Color::Modifier(Color::FG_RED), std::cerr);
LevelLogger debug("debug", Color::Modifier(Color::FG_MAGENTA));
LevelLogger timer("timer", Color::Modifier(Color::FG_DEFAULT));

MPI_Comm comm = MPI_COMM_WORLD;
std::string header_msg = "";


unsigned int & lvl() { return level; }
void set_comm(MPI_Comm new_comm)
{
	comm = new_comm;
}

void set_header_msg(std::string msg)
{
	header_msg = msg;
}

void init()
{
	log_level = std::make_unique<lmap_t>();
	std::vector<std::string> levels{"status", "info", "error", "memory", "debug", "timer"};

	int count = 0;
	for (auto lvl : levels) {
		(*log_level)[lvl] = 2<<count;
		count++;
	}

	std::vector<std::string> clevels = config::getvec<std::string>("log");

	for (auto clvl : clevels) level |= (*log_level)[clvl];

}


std::string header()
{
	std::ostringstream os;

	std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
	char buffer[256];
	Color::Modifier green(Color::FG_GREEN);
	Color::Modifier def(Color::FG_DEFAULT);

	strftime(buffer, sizeof(buffer), "%a %b %d %H:%M:%S %Y", &tm);
	os << green << "[BoxMG <" << buffer << ">" << header_msg << "] " << def;

	return os.str();
}

bool LevelLogger::active()
{
	if (!log_level) log::init();

	return (*log_level)[name] & level;
}

int LogLevelBuf::sync()
{
	if (!log_level) init();

	int mpi_init, mpi_final;
	MPI_Initialized(&mpi_init);
	MPI_Finalized(&mpi_final);
	int rank;
	if (!mpi_final and mpi_init)
		MPI_Comm_rank(comm, &rank);
	else rank = 0;

	if ((rank == 0) and (str() != "") and ((*log_level)[name] & level)) {
		Color::Modifier def(Color::FG_DEFAULT);
		ost << header() << color << str() << def;
		str("");
	} else str("");
	return !std::cout;
}

}}
