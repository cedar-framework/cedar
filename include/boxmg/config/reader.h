#ifndef CONFIG_READER_H
#define CONFIG_READER_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>

namespace boxmg {
namespace config
{
    class Reader
    {

        public:
            Reader();
            template <typename T>
	            Reader(T&& fname);
            template <typename OptionType>
                static void set(std::string path, OptionType val);
            template <typename OptionType>
                static void setvec(std::string path, std::vector<OptionType> vec);
            /**
             * Gets a config option by path.
             *
             * @param path Path to requested option in config file.
             * @return The requested option.
             */
            template <typename OptionType>
                static OptionType get(std::string path);

            template <typename OptionType>
                static OptionType get(std::string path, OptionType default_value);
            /**
             * Gets an array of config options by path.
             *
             * @param path Path to requested option in config file.
             * @return The requested options as a vector.
             */
            template <typename OptionType>
                static std::vector<OptionType> getvec(std::string path);
            /**
             * Sets a new config file for future Config::get calls.
             *
             * @param fname Location of new config file.
             */
            static void set_config(std::string fname);
            ~Reader() {instance = NULL;}

        private:
            void read();
            static Reader* instance;
            boost::property_tree::ptree pt;
            std::string fname;

    };

    template <typename T>
	    Reader::Reader(T&& fname): fname(std::forward<T>(fname)) { read(); }

    template <typename OptionType>
        void Reader::set(std::string path, OptionType val)
        {
            if (!instance)
                instance = new Reader();
            instance->pt.put(path, val);
        }

    template <typename OptionType>
        void Reader::setvec(std::string path, std::vector<OptionType> vec)
        {
            if (!instance)
                instance = new Reader();
            using boost::property_tree::ptree;
            ptree &pos = instance->pt.get_child(path);
            pos.clear();

            for (typename std::vector<OptionType>::iterator it=vec.begin();it != vec.end();++it){
                instance->pt.put(path + ".", *it);
                std::cout << *it << std::endl;
            }

        }

    template <typename OptionType>
        OptionType Reader::get(std::string path)
        {
            if (!instance)
                instance = new Reader();
            return instance->pt.get<OptionType>(path);
        }

    template <typename OptionType>
        OptionType Reader::get(std::string path, OptionType default_value)
        {
            if (!instance)
                instance = new Reader();
            return instance->pt.get<OptionType>(path, default_value);
        }

    template <typename OptionType>
        std::vector<OptionType> Reader::getvec(std::string path)
        {
            if (!instance)
                instance = new Reader();
            std::vector<OptionType> retvec;
            using boost::property_tree::ptree;
            boost::optional<ptree&> child = instance->pt.get_child_optional(path);
            if (child) {
                ptree &pos = instance->pt.get_child(path);
                std::for_each(pos.begin(), pos.end(), [&retvec](ptree::value_type &v) {
                        retvec.push_back(v.second.get_value<OptionType>());
                        });
            }
            return retvec;
        }

    template <typename OptionType>
        void set(std::string path, OptionType val)
        {
            Reader::set<OptionType>(path, val);
        }

    template <typename OptionType>
        void setvec(std::string path, std::vector<OptionType> vec)
        {
            Reader::setvec<OptionType>(path, vec);
        }

    template <typename OptionType>
        OptionType get(std::string path)
        {
            return Reader::get<OptionType>(path);
        }

    template <typename OptionType>
        OptionType get(std::string path, OptionType default_value)
        {
            return Reader::get<OptionType>(path, default_value);
        }

    template <typename OptionType>
        std::vector<OptionType> getvec(std::string path)
        {
            return Reader::getvec<OptionType>(path);
        }
}
}
#endif
