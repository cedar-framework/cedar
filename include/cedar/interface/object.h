#ifndef CEDAR_INTERFACE_OBJECT_H
#define CEDAR_INTERFACE_OBJECT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	int handle;
	int refcount;
	unsigned short kind;
	void *ptr;
} cedar_object;

#define CEDAR_KIND_MASK 0xf
#define CEDAR_KIND_NULL    0
#define CEDAR_KIND_TOPO    1
#define CEDAR_KIND_SOLVER  2
#define CEDAR_KIND_MAT     3
#define CEDAR_KIND_VEC     4

cedar_object *cedar_object_create(unsigned short kind);
int cedar_object_get(int handle, cedar_object **obj);
int cedar_object_incref(int handle);
int cedar_object_decref(int handle);
int cedar_object_free(int handle);


#ifdef __cplusplus
}
#endif

#endif
