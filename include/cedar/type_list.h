#ifndef CEDAR_TYPE_LIST_H
#define CEDAR_TYPE_LIST_H

#include <memory>
#include <map>
#include <tuple>

#include <cedar/util/log.h>

namespace cedar
{
	template<class... T> struct type_list {};

	template<class ta, class tb>
	struct type_cat;

	template<typename ... a, typename ... b>
	struct type_cat<type_list<a ...>, type_list<b ...>>
	{ typedef type_list<a ..., b ...> type; };


	template<class T1, class T2>
	struct same_type
	{
		static const bool value = false;
	};

	template<class T>
	struct same_type<T, T>
	{
		static const bool value = true;
	};


	template<class U>
	class type_map;

	template<typename... types>
	class type_map<type_list<types...>>
	{
	public:
		typedef std::tuple<std::map<std::string,
		                            std::shared_ptr<types>>...> mtype;
		mtype kerns;
		std::array<std::string, std::tuple_size<mtype>::value> defaults;

		template<int n, class T>
		struct map_of_type: same_type<std::shared_ptr<T>,
		                              typename std::tuple_element<n, mtype>::type::mapped_type>
		{};


		template<int n, class T, bool match = false>
		struct matching_index
		{
			static const std::size_t value = matching_index<n+1, T, map_of_type<n+1, T>::value>::value;
		};

		template<int n, class T>
		struct matching_index<n, T, true>
		{
			static const std::size_t value = n;
		};

		template<class T>
		std::shared_ptr<T> get_ptr()
		{
			const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
			if (i >= defaults.size())
				log::error << "<type_map> type not found" << std::endl;
			return std::get<i>(kerns)[defaults[i]];
		}


		template<class T>
		std::shared_ptr<T> get_ptr(const std::string & name)
		{
			const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
			if (i >= defaults.size())
				log::error << "<type_map> type not found" << std::endl;
			return std::get<i>(kerns)[name];
		}

		template<class T>
		T & get()
		{
			return *(get_ptr<T>());
		}


		template<class T>
		T & get(const std::string & name)
		{
			return *(get_ptr<T>(name));
		}


		template<class T, class rclass>
		void add(const std::string & name)
		{
			const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
			if (i >= defaults.size())
				log::error << "<type_map> type not found" << std::endl;
			std::get<i>(kerns)[name] = std::make_shared<rclass>();
		}


		template<class T, class rclass, class... Args>
		void add(const std::string & name, Args&&... args)
		{
			const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
			if (i >= defaults.size())
				log::error << "<type_map> type not found" << std::endl;
			std::get<i>(kerns)[name] = std::make_shared<rclass>(std::forward<Args>(args)...);
		}


		template<class T>
		void set(const std::string & name)
		{
			const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
			if (i >= defaults.size())
				log::error << "<type_map> type not found" << std::endl;
			defaults[i] = name;
		}
	};
}

#endif
