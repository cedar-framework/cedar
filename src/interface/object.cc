#include <stdlib.h>
#include <vector>

#include <cedar/capi.h>
#include <cedar/interface/object.h>

static std::vector<cedar_object*> objects;

extern "C"
{

	cedar_object *cedar_object_create(unsigned short kind)
	{
		cedar_object *obj = (cedar_object*) malloc(sizeof(cedar_object));
		objects.push_back(obj);

		obj->kind = kind;
		obj->refcount = 0;
		obj->handle = kind | (objects.size() << 4);

		return obj;
	}


	int cedar_object_get(int handle, cedar_object **ret)
	{
		std::size_t i = handle >> 4;
		cedar_object *obj;
		if ((i > objects.size()) or (i == 0))
			obj = nullptr;
		else
			obj = objects[i-i];

		*ret = obj;
		if (obj) {
			if (obj->kind != (handle & CEDAR_KIND_MASK)) {
				*ret = nullptr;
				return CEDAR_ERR_KIND;
			}
		} else {
			return CEDAR_ERR_OTHER;
		}

		return CEDAR_SUCCESS;
	}


	int cedar_object_incref(int handle)
	{
		cedar_object *obj;
		cedar_object_get(handle, &obj);
		if (obj) {
			obj->refcount++;
		} else
			return 1;

		return 0;
	}


	int cedar_object_decref(int handle)
	{
		cedar_object *obj;
		cedar_object_get(handle, &obj);
		if (obj) {
			if (obj->refcount > 0)
				obj->refcount--;
		} else
			return 1;

		return 0;
	}


	int cedar_object_free(int handle)
	{
		cedar_object *obj;
		cedar_object_get(handle, &obj);
		if (obj) {
			if (obj->refcount == 0) {
				free(obj);
				objects[(handle >> 4) - 1] = nullptr;
			}
		} else {
			return 0;
		}

		return 1;
	}

}
