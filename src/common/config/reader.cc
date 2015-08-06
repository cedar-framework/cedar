#include "reader.h"

#include <boost/filesystem.hpp>

namespace boxmg { namespace config
{
    Reader* Reader::instance = NULL;

    Reader::Reader():
        fname("config.json")
    {
        read();
    }

    void Reader::read()
    {
        if (!boost::filesystem::exists(fname))
            fname = "../" + fname;
        boost::property_tree::json_parser::read_json(fname, pt);
    }

    void Reader::set_config(std::string fname)
    {
	    if (!instance) {
            instance = new Reader(fname);
	    } else {
		    instance->fname = fname;
		    instance->read();
	    }
    }
}}
