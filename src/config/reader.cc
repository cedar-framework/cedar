#include <cedar/config/reader.h>

namespace cedar { namespace config
{
	reader::reader():
        fname("config.json")
    {
        read();
    }

    void reader::read()
    {
	    if (not (fname == "")) {
		    std::ifstream cfile(fname);
		    root = json::parse(cfile);
	    }
    }

    void reader::set_config(std::string fname)
    {
	    this->fname = fname;
    }

	std::shared_ptr<reader> reader::getconf(std::string path)
	{
	    std::replace(path.begin(), path.end(), '.', '/');
	    std::string rpath("/" + path);
	    try {
		    json stree = root[json::json_pointer(rpath)];
		    return std::make_shared<reader>(std::move(stree));
	    } catch(...) {
		    return nullptr;
	    }
	}
}}
