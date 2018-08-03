#include <cedar/config.h>

namespace cedar {

	config::config():
        fname("config.json")
    {
        read();
    }

    void config::read()
    {
	    if (not (fname == "")) {
		    std::ifstream cfile(fname);
		    root = json::parse(cfile);
	    }
    }

    void config::set_config(std::string fname)
    {
	    this->fname = fname;
    }

	std::shared_ptr<config> config::getconf(std::string path)
	{
	    std::replace(path.begin(), path.end(), '.', '/');
	    std::string rpath("/" + path);

	    json stree = (path == "") ? root : root[json::json_pointer(rpath)];

	    if ((stree.size() == 0) and (path != ""))
		    return nullptr;
	    else
		    return std::make_shared<config>(std::move(stree));
	}
}
