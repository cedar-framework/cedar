#include <memory>
#include <fstream>

#include <cedar/capi.h>
#include <cedar/interface/object.h>
#include <cedar/interface/config.h>


std::shared_ptr<cedar::config> cedar_config_getobj(cedar_config handle)
{
	if (CEDAR_GET_KIND(handle) != CEDAR_KIND_CONFIG)
		return nullptr;

	cedar_object *obj;
	int ierr = cedar_object_get(handle, &obj);
	if (ierr or (not obj))
		return nullptr;

	auto ptrptr = reinterpret_cast<std::shared_ptr<cedar::config>*>(obj->ptr);
	return *ptrptr;
}


static bool file_exists(const char *fname)
{
	std::ifstream tfile(fname);
	return tfile.good();
}


extern "C"
{
	int cedar_config_create(const char *fname, cedar_config *newconfig)
	{
		using namespace cedar;

		if (not file_exists(fname)) {
			*newconfig = CEDAR_CONFIG_NULL;
			return CEDAR_ERR_FNAME;
		}

		cedar_object *obj = cedar_object_create(CEDAR_KIND_CONFIG);

		std::shared_ptr<config> *conf_ptr;
		auto conf = std::make_shared<config>(fname);
		conf_ptr = new std::shared_ptr<config>(std::move(conf));
		obj->ptr = reinterpret_cast<void*>(conf_ptr);
		*newconfig = obj->handle;

		return CEDAR_SUCCESS;
	}


	int cedar_config_free(cedar_config *conf)
	{
		if (CEDAR_GET_KIND(*conf) != CEDAR_KIND_CONFIG)
			return CEDAR_ERR_CONFIG;

		cedar_object *obj;
		cedar_object_get(*conf, &obj);
		if (obj) {
			if (obj->refcount == 0) {
				auto *ptr = reinterpret_cast<std::shared_ptr<cedar::config>*>(obj->ptr);
				delete ptr;
				cedar_object_free(*conf);
			}
		} else
			return CEDAR_ERR_CONFIG;

		*conf = CEDAR_CONFIG_NULL;

		return CEDAR_SUCCESS;
	}
}
