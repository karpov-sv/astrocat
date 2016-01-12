#ifndef ASTROCAT_H
#define ASTROCAT_H

#ifndef MIN
#define MIN(a, b)				\
    ({ typeof (a) _a = (a);     \
        typeof (b) _b = (b);    \
        _a < _b ? _a : _b; })
#endif /* MIN */
#ifndef MAX
#define MAX(a, b) \
    ({ typeof (a) _a = (a);     \
        typeof (b) _b = (b);    \
        _a > _b ? _a : _b; })
#endif /* MAX */

typedef struct {
	float8 x;
	float8 y;
	float8 z;
} skypoint;

typedef struct {
	float8 x;
	float8 y;
	float8 z;
	float8 sr; /* SR is the linear distance! */
} skycircle;

typedef struct {
	float8 x[2];
	float8 y[2];
	float8 z[2];
} skybox;

Datum		skypoint_from_alpha_delta(PG_FUNCTION_ARGS);
Datum		skypoint_in(PG_FUNCTION_ARGS);
Datum		skypoint_out(PG_FUNCTION_ARGS);

Datum		skycircle_from_alpha_delta(PG_FUNCTION_ARGS);
Datum		skycircle_in(PG_FUNCTION_ARGS);
Datum		skycircle_out(PG_FUNCTION_ARGS);

Datum		skypoint_in_skycircle(PG_FUNCTION_ARGS);

Datum		skybox_in(PG_FUNCTION_ARGS);
Datum		skybox_out(PG_FUNCTION_ARGS);

skybox     *skybox_copy(skybox *);

Datum       skybox_compress(PG_FUNCTION_ARGS);
Datum       skybox_decompress(PG_FUNCTION_ARGS);
Datum       skybox_consistent(PG_FUNCTION_ARGS);
Datum       skybox_union(PG_FUNCTION_ARGS);
Datum       skybox_penalty(PG_FUNCTION_ARGS);
Datum       skybox_picksplit(PG_FUNCTION_ARGS);

#endif /* ASTROCAT_H */
