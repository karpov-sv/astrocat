CREATE EXTENSION astrocat;

-- test data type accuracy down to ~0.00036 arcsec
SELECT skypoint(1.0, 1.0);
SELECT skypoint(1.1, 1.1);
SELECT skypoint(1.01, 1.01);
SELECT skypoint(1.001, 1.001);
SELECT skypoint(1.0001, 1.0001);
SELECT skypoint(1.00001, 1.00001);
SELECT skypoint(1.000001, 1.000001);
SELECT skypoint(1.0000001, 1.0000001);

SELECT skycircle(0.0, 0.0, 0.1);
SELECT skycircle(0.0, 0.0, 1.0);
SELECT skycircle(0.0, 0.0, 1.1);
SELECT skycircle(0.0, 0.0, 1.01);
SELECT skycircle(0.0, 0.0, 1.001);
SELECT skycircle(0.0, 0.0, 1.0001);
SELECT skycircle(0.0, 0.0, 1.00001);
SELECT skycircle(0.0, 0.0, 1.000001);

-- test containment
SELECT skypoint(0.0, 0.0) @ skycircle(0.0, 0.0, 0.0);
SELECT skypoint(0.1, 0.0) @ skycircle(0.0, 0.0, 0.1);
SELECT skypoint(0.01, 0.0) @ skycircle(0.0, 0.0, 0.01);
SELECT skypoint(0.001, 0.0) @ skycircle(0.0, 0.0, 0.001);
SELECT skypoint(0.0001, 0.0) @ skycircle(0.0, 0.0, 0.0001);
SELECT skypoint(0.00001, 0.0) @ skycircle(0.0, 0.0, 0.00001);

SELECT skypoint(1.0, 1.0) @ skycircle(0.0, 0.0, 1.0);
SELECT skypoint(1.0, 1.0) @ skycircle(0.0, 0.0, 1.4);
SELECT skypoint(1.0, 1.0) @ skycircle(0.0, 0.0, 1.5);
SELECT skypoint(1.0, 1.0) @ skycircle(0.0, 0.0, 10.0);

-- Distances
SELECT *, skypoint(ra,dec) <-> skypoint(3,2) AS dist FROM (
SELECT generate_series(0,10,1) AS ra, generate_series(10,0,-1) AS dec) t
ORDER BY skypoint(ra,dec) <-> skypoint(3,2);
