#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <nlohmann/json.hpp>

namespace cedar {
namespace config
{
	using json = nlohmann::json;
    class reader
    {
        public:
	    reader();
    reader(json && stree) : root(std::move(stree)), fname("") {}
            template <typename T>
	            reader(T&& fname);
            template <typename OptionType>
                void set(std::string path, OptionType val);
            template <typename OptionType>
                void setvec(std::string path, std::vector<OptionType> vec);
            /**
             * Gets a config option by path.
             *
             * @param path Path to requested option in config file.
             * @return The requested option.
             */
            template <typename OptionType>
                OptionType get(std::string path);

            template <typename OptionType>
                OptionType get(std::string path, OptionType default_value);
            /**
             * Gets an array of config options by path.
             *
             * @param path Path to requested option in config file.
             * @return The requested options as a vector.
             */
            template <typename OptionType>
	            std::vector<OptionType> getvec(std::string path);
            template <typename OptionType>
	            std::vector<std::vector<OptionType>> getnvec(std::string path);
            std::shared_ptr<reader> getconf(std::string path);
            /**
             * Sets a new config file for future Config::get calls.
             *
             * @param fname Location of new config file.
             */
            void set_config(std::string fname);

        private:
            void read();
            json root;
            std::string fname;

    };

    template <typename T>
	    reader::reader(T&& fname): fname(std::forward<T>(fname)) { read(); }

    template <typename OptionType>
        void reader::set(std::string path, OptionType val)
        {
	        std::replace(path.begin(), path.end(), '.', '/');
	        std::string rpath("/" + path);
	        root[json::json_pointer(rpath)] = val;
        }

    template <typename OptionType>
        void reader::setvec(std::string path, std::vector<OptionType> vec)
        {
	        set<std::vector<OptionType>>(path, vec);
        }

    template <class otype>
        otype reader::get(std::string path)
        {
	        std::replace(path.begin(), path.end(), '.', '/');
	        std::string rpath("/" + path);
	        return root[json::json_pointer(rpath)].get<otype>();
        }
    template <class otype>
        otype reader::get(std::string path, otype default_value)
        {
	        std::replace(path.begin(), path.end(), '.', '/');
	        std::string rpath("/" + path);
	        try {
		        return root[json::json_pointer(rpath)].get<otype>();
	        } catch(...) {
		        return default_value;
	        }
        }

    template <class otype>
        std::vector<otype> reader::getvec(std::string path)
        {
	        try {
		        return get<std::vector<otype>>(path);
	        } catch(...) {
		        return std::vector<otype>();
	        }
        }

    template <typename otype>
	    std::vector<std::vector<otype>> reader::getnvec(std::string path)
    {
	    return get<std::vector<std::vector<otype>>>(path);
    }
}
}
#endif
