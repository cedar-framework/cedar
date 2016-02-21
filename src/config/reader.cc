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
        if (!boost::filesystem::exists(fname))
            fname = "../" + fname;
        boost::property_tree::json_parser::read_json(fname, pt);
    }

    void reader::set_config(std::string fname)
    {
	    this->fname = fname;
    }
}}
