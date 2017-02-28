#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>
#include <memory>

namespace boxmg {
namespace config
{
    class reader
    {
        public:
            reader();
            reader(boost::property_tree::ptree && ptri): pt(std::move(ptri)), fname("") {}
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
            boost::property_tree::ptree pt;
            std::string fname;

    };

    template <typename T>
	    reader::reader(T&& fname): fname(std::forward<T>(fname)) { read(); }

    template <typename OptionType>
        void reader::set(std::string path, OptionType val)
        {
	        pt.put(path, val);
        }

    template <typename OptionType>
        void reader::setvec(std::string path, std::vector<OptionType> vec)
        {
            using boost::property_tree::ptree;
            ptree &pos = pt.get_child(path);
            pos.clear();

            for (typename std::vector<OptionType>::iterator it=vec.begin();it != vec.end();++it){
                pt.put(path + ".", *it);
            }

        }

    template <typename OptionType>
        OptionType reader::get(std::string path)
        {
            return pt.get<OptionType>(path);
        }

    template <typename OptionType>
        OptionType reader::get(std::string path, OptionType default_value)
        {
            return pt.get<OptionType>(path, default_value);
        }

    template <typename OptionType>
        std::vector<OptionType> reader::getvec(std::string path)
        {
            std::vector<OptionType> retvec;
            using boost::property_tree::ptree;
            boost::optional<ptree&> child = pt.get_child_optional(path);
            if (child) {
                ptree &pos = pt.get_child(path);
                std::for_each(pos.begin(), pos.end(), [&retvec](ptree::value_type &v) {
                        retvec.push_back(v.second.get_value<OptionType>());
                        });
            }
            return retvec;
        }

    template <typename OptionType>
	    std::vector<std::vector<OptionType>> reader::getnvec(std::string path)
    {
	    std::vector<std::vector<OptionType>> retvec;
	    using boost::property_tree::ptree;
	    boost::optional<ptree&> child = pt.get_child_optional(path);
	    if (child) {
		    ptree &pos = pt.get_child(path);
		    std::for_each(pos.begin(), pos.end(), [&retvec](ptree::value_type &v) {
				    std::vector<OptionType> toadd;
				    std::for_each(v.second.begin(), v.second.end(), [&toadd](ptree::value_type &v) {
						    toadd.push_back(v.second.get_value<OptionType>());
					    });
				    retvec.push_back(toadd);
			    });
	    }
	    return retvec;
    }
}
}
#endif
