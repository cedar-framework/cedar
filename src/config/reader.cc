#include <boxmg/config/reader.h>

#include <boost/filesystem.hpp>

namespace boxmg { namespace config
{
	reader::reader():
        fname("config.json")
    {
        read();
    }

    void reader::read()
    {
	    if (!(fname == "")) {
		    if (!boost::filesystem::exists(fname))
			    fname = "../" + fname;
		    boost::property_tree::json_parser::read_json(fname, pt);
	    }
    }

    void reader::set_config(std::string fname)
    {
	    this->fname = fname;
    }

	std::shared_ptr<reader> reader::getconf(std::string path)
	{
		using boost::property_tree::ptree;

		boost::optional<ptree&> child = pt.get_child_optional(path);
		if (child) {
			ptree kid = pt.get_child(path);
			return std::make_shared<reader>(std::move(kid));
		} else {
			return nullptr;
		}
	}
}}
