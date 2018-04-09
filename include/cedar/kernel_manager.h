#ifndef CEDAR_KERNEL_MANAGER_H
#define CEDAR_KERNEL_MANAGER_H

#include <memory>
#include <map>
#include <tuple>

#include <cedar/type_list.h>
#include <cedar/kernel_params.h>
#include <cedar/config/reader.h>


namespace cedar {

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
class kernel_manager;

template<typename... types>
class kernel_manager<type_list<types...>>
{
public:
    kernel_manager(std::shared_ptr<kernel_params> params) : params(params) {}
	kernel_manager(config::reader & conf) { params = build_kernel_params(conf); }

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
	T & get()
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		return *std::get<i>(kerns)[defaults[i]];
	}


	template<class T>
	T & get(const std::string & name)
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		return *std::get<i>(kerns)[name];
	}


	template<class T, class rclass>
	void add(const std::string & name)
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		std::get<i>(kerns)[name] = std::make_shared<rclass>();
		this->init<T>(name);
	}


	template<class T, class rclass, class... Args>
	void add(const std::string & name, Args&&... args)
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		std::get<i>(kerns)[name] = std::make_shared<rclass>(std::forward<Args>(args)...);
		this->init<T>(name);
	}


	template<class T>
	void set(const std::string & name)
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		defaults[i] = name;
	}


	template<class T, class... Args>
	void setup(Args&&... args)
	{
		this->get<T>().setup(std::forward<Args>(args)...);
	}

	template<class T, class... Args>
	void run(Args&&... args)
	{
		const std::size_t i = matching_index<0, T, map_of_type<0, T>::value>::value;
		std::get<i>(kerns)[defaults[i]]->run(std::forward<Args>(args)...);
	}


	template<class T>
	void init(const std::string & name)
	{
		this->get<T>(name).add_params(params);
	}

protected:
	std::shared_ptr<kernel_params> params;
};

}

#endif
