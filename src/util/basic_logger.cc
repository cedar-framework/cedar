#include <cedar/util/basic_logger.h>

using namespace cedar;

basic_logger::basic_logger(global_params & params)
	: level(params.log_level), header_msg(""), comm(MPI_COMM_WORLD),
	  status(loglevel::status, level, header_msg, comm, Color::Modifier(Color::FG_DEFAULT)),
	  error(loglevel::error, level, header_msg, comm, Color::Modifier(Color::FG_RED), std::cerr),
	  info(loglevel::info, level, header_msg, comm, Color::Modifier(Color::FG_DEFAULT)),
	  debug(loglevel::debug, level, header_msg, comm, Color::Modifier(Color::FG_MAGENTA)),
	  timer(loglevel::timer, level, header_msg, comm, Color::Modifier(Color::FG_DEFAULT)),
	  memory(loglevel::memory, level, header_msg, comm, Color::Modifier(Color::FG_BLUE)) {}



void basic_logger::push_comm(MPI_Comm new_comm)
{
	saved_comms.push(comm);
	comm = new_comm;
}


void basic_logger::pop_comm()
{
	comm = saved_comms.top();
	saved_comms.pop();
}


void basic_logger::push_level(std::string header, loglevel_t lvl)
{
	saved_levels.emplace(level, header_msg);
	level = lvl;
	set_header(" (" + header + " " + std::to_string(saved_levels.size()-1) + ")");
}


void basic_logger::pop_level()
{
	auto prev = saved_levels.top();
	level = std::get<0>(prev);
	set_header(std::get<1>(prev));
	saved_levels.pop();
}


std::string LogLevelBuf::header()
{
	std::ostringstream os;

	std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
	char buffer[256];
	Color::Modifier green(Color::FG_GREEN);
	Color::Modifier def(Color::FG_DEFAULT);

	strftime(buffer, sizeof(buffer), "%a %b %d %H:%M:%S %Y", &tm);
	os << green << "[Cedar <" << buffer << ">" << header_msg << "] " << def;

	return os.str();
}


int LogLevelBuf::sync()
{
	int mpi_init, mpi_final;
	MPI_Initialized(&mpi_init);
	MPI_Finalized(&mpi_final);
	int rank;
	if (!mpi_final and mpi_init)
		MPI_Comm_rank(comm, &rank);
	else rank = 0;

	if ((rank == 0) and (str() != "") and (curlevel & levelid)) {
		Color::Modifier def(Color::FG_DEFAULT);
		ost << header() << color << str() << def;
		str("");
	} else str("");
	return !std::cout;
}
