-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION astrocat" to load this file. \quit

-- Spherical point, I/O is in degrees
CREATE FUNCTION skypoint(FLOAT8, FLOAT8)
   RETURNS skypoint
   AS 'MODULE_PATHNAME' , 'skypoint_from_alpha_delta'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE FUNCTION skypoint_in(CSTRING)
   RETURNS skypoint
   AS 'MODULE_PATHNAME' , 'skypoint_in'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE FUNCTION skypoint_out(skypoint)
   RETURNS CSTRING
   AS 'MODULE_PATHNAME' , 'skypoint_out'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE TYPE skypoint (
   internallength = 24,
   input = skypoint_in,
   output = skypoint_out
);

CREATE FUNCTION skypoint_equal(skypoint,skypoint)
   RETURNS BOOL
   AS 'MODULE_PATHNAME' , 'skypoint_equal'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE OPERATOR  = (
   LEFTARG    = skypoint,
   RIGHTARG   = skypoint,
   COMMUTATOR = =,
   NEGATOR    = <>,
   PROCEDURE  = skypoint_equal,
   RESTRICT   = contsel,
   JOIN       = contjoinsel
);

CREATE FUNCTION skypoint_equal_neg(skypoint,skypoint)
   RETURNS BOOL
   AS 'MODULE_PATHNAME' , 'skypoint_equal_neg'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE OPERATOR  <> (
   LEFTARG    = skypoint,
   RIGHTARG   = skypoint,
   COMMUTATOR = <>,
   NEGATOR    = =,
   PROCEDURE  = skypoint_equal_neg ,
   RESTRICT   = contsel,
   JOIN       = contjoinsel
);

CREATE FUNCTION skypoint_distance(skypoint,skypoint)
   RETURNS FLOAT8
   AS 'MODULE_PATHNAME' , 'skypoint_distance'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE OPERATOR  <-> (
   LEFTARG    = skypoint,
   RIGHTARG   = skypoint,
   COMMUTATOR = '<->',
   PROCEDURE  = skypoint_distance
);

-- Spherical circle, I/O is in degrees

CREATE FUNCTION skycircle(FLOAT8, FLOAT8, FLOAT8)
   RETURNS skycircle
   AS 'MODULE_PATHNAME' , 'skycircle_from_alpha_delta'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE FUNCTION skycircle_in(CSTRING)
   RETURNS skycircle
   AS 'MODULE_PATHNAME' , 'skycircle_in'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE FUNCTION skycircle_out(skycircle)
   RETURNS CSTRING
   AS 'MODULE_PATHNAME' , 'skycircle_out'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE TYPE skycircle (
   internallength = 32,
   input = skycircle_in,
   output = skycircle_out
);

CREATE OR REPLACE FUNCTION skypoint_in_skycircle(skypoint, skycircle)
RETURNS BOOL
AS 'MODULE_PATHNAME'
LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION skypoint_in_skycircle(float8, float8, float8, float8, float8)
RETURNS BOOL
AS 'SELECT skypoint_in_skycircle(skypoint($1,$2),skycircle($3,$4,$5));'
LANGUAGE 'sql'
IMMUTABLE STRICT ;

CREATE OPERATOR @ (
   LEFTARG    = skypoint,
   RIGHTARG   = skycircle,
   PROCEDURE  = skypoint_in_skycircle,
   RESTRICT   = contsel,
   JOIN       = contjoinsel
);

-- 3D box, to be internally used as a GiST key

CREATE FUNCTION skybox_in(CSTRING)
   RETURNS skybox
   AS 'MODULE_PATHNAME' , 'skybox_in'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE FUNCTION skybox_out(skybox)
   RETURNS CSTRING
   AS 'MODULE_PATHNAME' , 'skybox_out'
   LANGUAGE 'c'
   IMMUTABLE STRICT ;

CREATE TYPE skybox (
   internallength = 48,
   input = skybox_in,
   output = skybox_out
);

CREATE OR REPLACE FUNCTION skybox_compress(internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_decompress(internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_consistent(internal, internal, smallint, oid, internal)
RETURNS bool
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_union(internal, internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_penalty(internal, internal, internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_picksplit(internal, internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE OR REPLACE FUNCTION skybox_same(internal, internal, internal)
RETURNS internal
AS 'MODULE_PATHNAME'
LANGUAGE C STRICT;

CREATE FUNCTION skybox_distance(internal, internal, int4, oid)
RETURNS float8
AS 'MODULE_PATHNAME'
LANGUAGE 'c';

-- GiST operator class

CREATE OPERATOR CLASS skybox
   DEFAULT FOR TYPE skypoint USING gist AS
--   OPERATOR   1 = (skypoint, skypoint),
   OPERATOR  11 @ (skypoint, skycircle),
   OPERATOR  100 <-> (skypoint, skypoint) FOR ORDER BY float_ops,
   FUNCTION  1 skybox_consistent (internal, internal, smallint, oid, internal),
   FUNCTION  2 skybox_union (internal, internal),
   FUNCTION  3 skybox_compress (internal),
   FUNCTION  4 skybox_decompress (internal),
   FUNCTION  5 skybox_penalty (internal, internal, internal),
   FUNCTION  6 skybox_picksplit (internal, internal),
   FUNCTION  7 skybox_same (internal, internal, internal),
   FUNCTION  8 skybox_distance (internal, internal, int4, oid),
   STORAGE   skybox;
