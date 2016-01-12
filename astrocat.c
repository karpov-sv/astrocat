/*-------------------------------------------------------------------------
 *
 * astrocat.c
 *     astro catalogue extenstion for PostgreSQL
 *
 * Copyright (c) 2014-2015, PostgreSQL Global Development Group
 *
 * IDENTIFICATION
 *		hello_ext/hello_ext.c
 *
 *-------------------------------------------------------------------------
 */

#include <math.h>

#include "postgres.h"
#include "fmgr.h"
#include "utils/builtins.h"
#include "access/gist.h"
#include "access/skey.h"

#include "astrocat.h"

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(skypoint_from_alpha_delta);
PG_FUNCTION_INFO_V1(skypoint_in);
PG_FUNCTION_INFO_V1(skypoint_out);
PG_FUNCTION_INFO_V1(skypoint_equal);
PG_FUNCTION_INFO_V1(skypoint_equal_neg);
PG_FUNCTION_INFO_V1(skypoint_distance);

PG_FUNCTION_INFO_V1(skycircle_from_alpha_delta);
PG_FUNCTION_INFO_V1(skycircle_in);
PG_FUNCTION_INFO_V1(skycircle_out);

PG_FUNCTION_INFO_V1(skypoint_in_skycircle);

PG_FUNCTION_INFO_V1(skybox_in);
PG_FUNCTION_INFO_V1(skybox_out);

PG_FUNCTION_INFO_V1(skybox_compress);
PG_FUNCTION_INFO_V1(skybox_decompress);
PG_FUNCTION_INFO_V1(skybox_union);
PG_FUNCTION_INFO_V1(skybox_same);
PG_FUNCTION_INFO_V1(skybox_consistent);
PG_FUNCTION_INFO_V1(skybox_penalty);
PG_FUNCTION_INFO_V1(skybox_picksplit);
PG_FUNCTION_INFO_V1(skybox_distance);

float8 distance_from_linear(float8 );
float8 linear_from_distance(float8 );

skybox *skybox_union_two(skybox *, skybox *);
skybox *skybox_inter_two (skybox *, skybox *);
double skybox_union_size(skybox *, skybox *);
double skybox_size(skybox *);

#define Sqr(a) ((a)*(a))

Datum
skypoint_from_alpha_delta(PG_FUNCTION_ARGS)
{
	skypoint *p = (skypoint *) palloc(sizeof(skypoint));
	float8 ra = PG_GETARG_FLOAT8(0);
	float8 dec = PG_GETARG_FLOAT8(1);

	p->x = cos(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->y = sin(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->z = sin(dec*M_PI/180.0);

	PG_RETURN_POINTER(p);
}

Datum
skypoint_in(PG_FUNCTION_ARGS)
{
	skypoint *p = (skypoint *) palloc(sizeof(skypoint));
	char *c = PG_GETARG_CSTRING(0);
	float8 ra = 0;
	float8 dec = 0;

	sscanf(c, "(%lf,%lf)", &ra, &dec);

	p->x = cos(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->y = sin(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->z = sin(dec*M_PI/180.0);

	PG_RETURN_POINTER(p);
}

Datum
skypoint_out(PG_FUNCTION_ARGS)
{
	skypoint *p = (skypoint *) PG_GETARG_POINTER(0);
	char *buffer = (char *) palloc(1024);
	float8 ra = 0;
	float8 dec = 0;

	dec = asin(p->z)*180.0/M_PI;
	ra = atan2(p->y, p->x)*180.0/M_PI;

	while(ra < 0)
		ra += 360.0;

	sprintf(buffer, "(%.9f,%.9f)", ra, dec);

	PG_RETURN_CSTRING(buffer);
}

Datum
skypoint_equal(PG_FUNCTION_ARGS)
{
	skypoint *p1 = (skypoint *) PG_GETARG_POINTER(0);
	skypoint *p2 = (skypoint *) PG_GETARG_POINTER(1);

	PG_RETURN_BOOL( memcmp(p1, p2, sizeof(skypoint)) ? false : true);
}

Datum
skypoint_equal_neg(PG_FUNCTION_ARGS)
{
	skypoint *p1 = (skypoint *) PG_GETARG_POINTER(0);
	skypoint *p2 = (skypoint *) PG_GETARG_POINTER(1);

	PG_RETURN_BOOL( !memcmp(p1, p2, sizeof(skypoint)) ? false : true );
}

inline float8 distance_from_linear(float8 dist)
{
#ifdef UNNECESSARY_OPTIMIZATION
	if (dist < 0.5)
		/* Approximate conversion of linear distance to angular one */
		/* Accuracy is ~1 arcsec at 30 deg, ~170 arcsec at 60 deg and 2300 arcsec at 90 deg */
		return 180.0*dist/M_PI + 15.0*dist*dist*dist/(2*M_PI) + 27.0*dist*dist*dist*dist*dist/(32*M_PI);
	else
#endif
		/* Exact formula is: 180*arctan(1/2*(-x^2+4)^(1/2)*x,-1/2*x^2+1)/Pi */
		return 180.0*atan2(0.5*sqrt(4.0 - dist*dist)*dist, -0.5*dist*dist + 1)/M_PI;
}

inline float8 linear_from_distance(float8 sr)
{
#ifdef UNNECESSARY_OPTIMIZATION
	if (sr < 30)
		return M_PI/180.0*sr - M_PI*M_PI*M_PI/139968000.0*sr*sr*sr + M_PI*M_PI*M_PI*M_PI*M_PI/362797056000000.0*sr*sr*sr*sr*sr;
	else
#endif
		return sqrt(Sqr(1.0 - cos(sr*M_PI/180.0)) + Sqr(sin(sr*M_PI/180.0)));
}

inline float8 skypoint_linear_dist_squared(skypoint *p1, skypoint *p2)
{
	float8 dx = p1->x - p2->x;
	float8 dy = p1->y - p2->y;
	float8 dz = p1->z - p2->z;

	return dx*dx + dy*dy + dz*dz;
}

inline float8 skypoint_dist(skypoint *p1, skypoint *p2)
{
	float8 dist = sqrt(skypoint_linear_dist_squared(p1, p2)); /* Linear distance */

	return distance_from_linear(dist);
}

Datum
skypoint_distance(PG_FUNCTION_ARGS)
{
	skypoint *p1 = (skypoint *) PG_GETARG_POINTER(0);
	skypoint *p2 = (skypoint *) PG_GETARG_POINTER(1);

	PG_RETURN_FLOAT8( skypoint_dist(p1, p2) );
}

Datum
skycircle_from_alpha_delta(PG_FUNCTION_ARGS)
{
	skycircle *p = (skycircle *) palloc(sizeof(skycircle));
	float8 ra = PG_GETARG_FLOAT8(0);
	float8 dec = PG_GETARG_FLOAT8(1);
	float8 sr = PG_GETARG_FLOAT8(2);

	p->x = cos(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->y = sin(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->z = sin(dec*M_PI/180.0);
	p->sr = linear_from_distance(sr);

	PG_RETURN_POINTER(p);
}

Datum
skycircle_in(PG_FUNCTION_ARGS)
{
	skycircle *p = (skycircle *) palloc(sizeof(skycircle));
	char *c = PG_GETARG_CSTRING(0);
	float8 ra = 0;
	float8 dec = 0;
	float8 sr = 0;

	sscanf(c, "(%lf,%lf,%lf)", &ra, &dec, &sr);

	p->x = cos(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->y = sin(ra*M_PI/180.0)*cos(dec*M_PI/180.0);
	p->z = sin(dec*M_PI/180.0);
	p->sr = linear_from_distance(sr);

	PG_RETURN_POINTER(p);
}

Datum
skycircle_out(PG_FUNCTION_ARGS)
{
	skycircle *p = (skycircle *) PG_GETARG_POINTER(0);
	char *buffer = (char *) palloc(1024);
	float8 ra = 0;
	float8 dec = 0;
	float8 sr = 0;

	dec = asin(p->z)*180.0/M_PI;
	ra = atan2(p->y, p->x)*180.0/M_PI;
	sr = distance_from_linear(p->sr);

	while(ra < 0)
		ra += 360.0;

	sprintf(buffer, "(%.9f,%.9f,%.9f)", ra, dec, sr);

	PG_RETURN_CSTRING(buffer);
}

Datum
skypoint_in_skycircle(PG_FUNCTION_ARGS)
{
	skypoint *p = (skypoint *) PG_GETARG_POINTER(0);
	skycircle *c = (skycircle *) PG_GETARG_POINTER(1);
	bool result = false;

	if ( skypoint_linear_dist_squared(p, (skypoint *)c) <= c->sr*c->sr )
		result = true;

	PG_RETURN_BOOL(result);
}

Datum
skybox_in(PG_FUNCTION_ARGS)
{
	skybox *p = (skybox *) palloc(sizeof(skybox));
	char *c = PG_GETARG_CSTRING(0);

	sscanf(c, "(%lf,%lf,%lf),(%lf,%lf,%lf)",
		   &p->x[0], &p->y[0], &p->z[0],
		   &p->x[1], &p->y[1], &p->z[1]);

	PG_RETURN_POINTER(p);
}

Datum
skybox_out(PG_FUNCTION_ARGS)
{
	skybox *p = (skybox *) PG_GETARG_POINTER(0);
	char *buffer = (char *) palloc(1024);

	sprintf(buffer, "(%.9f,%.9f,%.9f),(%.9f,%.9f,%.9f)",
			p->x[0], p->y[0], p->z[0],
			p->x[1], p->y[1], p->z[1]);

	PG_RETURN_CSTRING(buffer);
}

Datum
skybox_compress(PG_FUNCTION_ARGS) {
	GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	GISTENTRY  *retval;

	if (entry->leafkey)
	{
		retval = palloc(sizeof(GISTENTRY));

		if (DatumGetPointer(entry->key) != NULL)
		{
			skybox *key = (skybox *) palloc(sizeof(skybox));
			skypoint *p = (skypoint *) DatumGetPointer(entry->key);

			key->x[0] = p->x;
			key->y[0] = p->y;
			key->z[0] = p->z;

			key->x[1] = p->x;
			key->y[1] = p->y;
			key->z[1] = p->z;

			gistentryinit(*retval, PointerGetDatum(key), entry->rel,
						  entry->page, entry->offset, false);
		}
		else
		{
			gistentryinit(*retval, (Datum) 0, entry->rel, entry->page,
						  entry->offset, false);
		}
	}
	else
	{
		retval = entry;
	}

	PG_RETURN_POINTER(retval);
}

Datum
skybox_decompress(PG_FUNCTION_ARGS)
{
	PG_RETURN_DATUM(PG_GETARG_DATUM(0));
}

Datum
skybox_same(PG_FUNCTION_ARGS)
{
	skybox *p1 = (skybox *)PG_GETARG_POINTER(0);
	skybox *p2 = (skybox *)PG_GETARG_POINTER(1);
	bool	   *result = (bool *) PG_GETARG_POINTER(2);

	*result = memcmp(p1, p2, sizeof(skybox)) ? false : true;

	PG_RETURN_POINTER(result);
}

Datum
skybox_consistent(PG_FUNCTION_ARGS)
{
    GISTENTRY  *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
    void  *query = PG_GETARG_POINTER(1);
    StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2);
    /* Oid subtype = PG_GETARG_OID(3); */
    bool *recheck = (bool *) PG_GETARG_POINTER(4);
    skybox  *key = (skybox *) DatumGetPointer(entry->key);
    bool retval = false;

    *recheck = false;        /* or false if check is exact */

	if (!GIST_LEAF(entry))
	{
		if (strategy == 11)
		{
			skycircle *c = (skycircle *)query;

			if(c->x >= key->x[0] - c->sr && c->x <= key->x[1] + c->sr &&
			   c->y >= key->y[0] - c->sr && c->y <= key->y[1] + c->sr &&
			   c->z >= key->z[0] - c->sr && c->z <= key->z[1] + c->sr)
				retval = true;
		}
	}
	else
	{
		if (strategy == 11)
		{
			skycircle *c = (skycircle *)query;
			float8 dx = key->x[0] - c->x;
			float8 dy = key->y[0] - c->y;
			float8 dz = key->z[0] - c->z;

			/* We use linear distance here to speed things up */
			if(dx*dx + dy*dy + dz*dz < c->sr*c->sr)
				retval = true;
		}
	}

    PG_RETURN_BOOL(retval);
}

skybox *skybox_copy(skybox *in)
{
	skybox *out = (skybox *) palloc(sizeof(skybox));

	memcpy(out, in, sizeof(skybox));

	return out;
}

Datum
skybox_union(PG_FUNCTION_ARGS)
{
    GistEntryVector *entryvec = (GistEntryVector *) PG_GETARG_POINTER(0);
    GISTENTRY *ent = entryvec->vector;
    skybox *out;
    int numranges;
	int i = 0;

    numranges = entryvec->n;
	out = skybox_copy((skybox *)DatumGetPointer(ent[0].key));

    for (i = 1; i < numranges; i++)
    {
		skybox *tmp = (skybox *) DatumGetPointer(ent[i].key);

		skybox_union_two(out, tmp);
    }

    PG_RETURN_POINTER(out);
}

Datum
skybox_penalty(PG_FUNCTION_ARGS)
{
    GISTENTRY *origentry = (GISTENTRY *) PG_GETARG_POINTER(0);
    GISTENTRY *newentry = (GISTENTRY *) PG_GETARG_POINTER(1);
    float *penalty = (float *) PG_GETARG_POINTER(2);
    skybox *orig = (skybox *) DatumGetPointer(origentry->key);
    skybox *new = (skybox *) DatumGetPointer(newentry->key);

    *penalty = skybox_union_size(orig, new) - skybox_size(orig);

    PG_RETURN_POINTER(penalty);
}

skybox *
skybox_union_two(skybox *a, skybox *b)
{
	a->x[0] = MIN(a->x[0], b->x[0]);
	a->y[0] = MIN(a->y[0], b->y[0]);
	a->z[0] = MIN(a->z[0], b->z[0]);

	a->x[1] = MAX(a->x[1], b->x[1]);
	a->y[1] = MAX(a->y[1], b->y[1]);
	a->z[1] = MAX(a->z[1], b->z[1]);

	return a;
}

double
skybox_union_size(skybox *a, skybox *b)
{
	float8 dx = MAX(a->x[1], b->x[1]) - MIN(a->x[0], b->x[0]);
	float8 dy = MAX(a->y[1], b->y[1]) - MIN(a->y[0], b->y[0]);
	float8 dz = MAX(a->z[1], b->z[1]) - MIN(a->z[0], b->z[0]);

	return dx*dy*dz;
}

double
skybox_size(skybox *a)
{
	float8 dx = a->x[1] - a->x[0];
	float8 dy = a->y[1] - a->y[0];
	float8 dz = a->z[1] - a->z[0];

	return dx*dy*dz;
}

skybox *
skybox_inter_two (skybox *a , skybox *b ) {
	if (a->x[1] < b->x[0] || b->x[1] < a->x[0]) return NULL;
	if (a->y[1] < b->y[0] || b->y[1] < a->y[0]) return NULL;
	if (a->z[1] < b->z[0] || b->z[1] < a->z[0]) return NULL;

	a->x[0] = MAX(a->x[0], b->x[0]);
	a->y[0] = MAX(a->y[0], b->y[0]);
	a->z[0] = MAX(a->z[0], b->z[0]);

	a->x[1] = MIN(a->x[1], b->x[1]);
	a->y[1] = MIN(a->y[1], b->y[1]);
	a->z[1] = MIN(a->z[1], b->z[1]);

	return a;
}

/* The following picksplit code is copied from pgSphere 1.1.1 */

#define WISH_F(a,b,c) (double)( -(double)(((a)-(b))*((a)-(b))*((a)-(b)))*(c) )

typedef struct
{
	double      cost;
	OffsetNumber pos;
} SPLITCOST;

static int
comparecost(const void *a, const void *b)
{
    if (((SPLITCOST *) a)->cost == ((SPLITCOST *) b)->cost)
		return 0;
	else
		return (((SPLITCOST *) a)->cost > ((SPLITCOST *) b)->cost) ? 1 : -1;
}

Datum
skybox_picksplit(PG_FUNCTION_ARGS)
{
	GistEntryVector    *entryvec = ( GistEntryVector *) PG_GETARG_POINTER(0);
    GIST_SPLITVEC  *v = (GIST_SPLITVEC *) PG_GETARG_POINTER(1);
    OffsetNumber  i,j;
	skybox *datum_alpha, *datum_beta;
    skybox datum_l, datum_r;
    skybox union_dl,union_dr;
    skybox union_d, inter_d;
    double size_alpha, size_beta;
    double size_waste, waste=-1.0;
    double size_l, size_r;
    int nbytes;
    OffsetNumber	seed_1 = 0, seed_2 = 0;
    OffsetNumber *left, *right;
    OffsetNumber maxoff;
    SPLITCOST  *costvector;

	maxoff  = entryvec->n - 1;
    nbytes = (maxoff + 2) * sizeof(OffsetNumber);
    v->spl_left  = (OffsetNumber *) palloc(nbytes);
    v->spl_right = (OffsetNumber *) palloc(nbytes);

    for (i = FirstOffsetNumber; i < maxoff; i = OffsetNumberNext(i)) {
	 	datum_alpha = (skybox *) DatumGetPointer(entryvec->vector[i].key);

		for (j = OffsetNumberNext(i); j <= maxoff; j = OffsetNumberNext(j)) {
		 	datum_beta = (skybox *) DatumGetPointer(entryvec->vector[j].key);

			memcpy( (void*)&union_d, (void*)datum_alpha, sizeof(skybox) );
			memcpy( (void*)&inter_d, (void*)datum_alpha, sizeof(skybox) );

			size_waste = skybox_size( skybox_union_two(&union_d, datum_beta) ) -
				( (skybox_inter_two( &inter_d, datum_beta )) ? skybox_size(&inter_d) : 0 );

			if (size_waste > waste) {
				waste = size_waste;
				seed_1 = i;
				seed_2 = j;
			}
		}
	}

	left = v->spl_left;
	v->spl_nleft = 0;
	right = v->spl_right;
	v->spl_nright = 0;

	if (seed_1 == 0 || seed_2 == 0) {
		seed_1 = 1;
		seed_2 = 2;
	}

	memcpy( (void*)&datum_l, (void *) DatumGetPointer(entryvec->vector[seed_1].key), sizeof(skybox) );
	memcpy( (void*)&datum_r, (void *) DatumGetPointer(entryvec->vector[seed_2].key), sizeof(skybox) );
	size_l = skybox_size(&datum_l);
	size_r = skybox_size(&datum_r);

	costvector = (SPLITCOST *) palloc(sizeof(SPLITCOST) * maxoff);
	for (i = FirstOffsetNumber; i <= maxoff; i = OffsetNumberNext(i)) {
		costvector[i - 1].pos = i;
	 	datum_alpha = (skybox *) DatumGetPointer(entryvec->vector[i].key);
		memcpy( (void*)&union_dl, (void*)&datum_l, sizeof(skybox) );
		skybox_union_two(&union_dl, datum_alpha );
		memcpy( (void*)&union_dr, (void*)&datum_r, sizeof(skybox) );
		skybox_union_two(&union_dr, datum_alpha );
		costvector[i - 1].cost = fabs( (skybox_size(&union_dl) - size_l) - (skybox_size(&union_dr) - size_r) );
	}
	qsort((void *) costvector, maxoff, sizeof(SPLITCOST), comparecost);

	for (j = 0; j < maxoff; j++) {
		i = costvector[j].pos;

		if (i == seed_1) {
			*left++ = i;
			v->spl_nleft++;
			continue;
		} else if (i == seed_2) {
			*right++ = i;
			v->spl_nright++;
			continue;
		}

	 	datum_alpha = (skybox *) DatumGetPointer(entryvec->vector[i].key);
		memcpy( (void*)&union_dl, (void*)&datum_l, sizeof(skybox) );
		memcpy( (void*)&union_dr, (void*)&datum_r, sizeof(skybox) );

		skybox_union_two(&union_dl, datum_alpha );
		skybox_union_two(&union_dr, datum_alpha );

		size_alpha = skybox_size( &union_dl );
		size_beta  = skybox_size( &union_dr );

		/* pick which page to add it to */
		if (size_alpha - size_l < size_beta - size_r + WISH_F(v->spl_nleft, v->spl_nright, 1.0e-9) ) {
			memcpy( (void*)&datum_l, (void*)&union_dl, sizeof(skybox) );
			size_l = size_alpha;
			*left++ = i;
			v->spl_nleft++;
		} else {
			memcpy( (void*)&datum_r, (void*)&union_dr, sizeof(skybox) );
			size_r = size_beta;
			*right++ = i;
			v->spl_nright++;
		}
	}

	pfree(costvector);

	*right = *left = FirstOffsetNumber;

	datum_alpha = (skybox *)palloc(sizeof(skybox));
	memcpy( (void*)datum_alpha, (void*)&datum_l, sizeof(skybox) );
	datum_beta  = (skybox *)palloc(sizeof(skybox));
	memcpy( (void*)datum_beta, (void*)&datum_r, sizeof(skybox) );

	v->spl_ldatum = PointerGetDatum(datum_alpha);
	v->spl_rdatum = PointerGetDatum(datum_beta);

	PG_RETURN_POINTER(v);
}

/* End of picksplit code */

Datum
skybox_distance(PG_FUNCTION_ARGS)
{
	GISTENTRY *entry = (GISTENTRY *) PG_GETARG_POINTER(0);
	skypoint *query = (skypoint *) PG_GETARG_POINTER(1);

	/* StrategyNumber strategy = (StrategyNumber) PG_GETARG_UINT16(2); */
	skybox *key = (skybox *) DatumGetPointer(entry->key);

	if (GIST_LEAF(entry))
	{
		skypoint point;

		point.x = key->x[0];
		point.y = key->y[0];
		point.z = key->z[0];

		PG_RETURN_FLOAT8( skypoint_dist(&point, query) );
	}
	else
	{
		float8 sum = 0.0;

		if (query->x < key->x[0])
			sum += Sqr(query->x - key->x[0]);
		else if (query->x > key->x[1])
			sum += Sqr(query->x - key->x[1]);

		if (query->y < key->y[0])
			sum += Sqr(query->y - key->y[0]);
		else if (query->y > key->x[1])
			sum += Sqr(query->y - key->y[1]);

		if (query->z < key->z[0])
			sum += Sqr(query->z - key->z[0]);
		else if (query->z > key->z[1])
			sum += Sqr(query->z - key->z[1]);

		PG_RETURN_FLOAT8( sqrt(sum) );
	}
}
