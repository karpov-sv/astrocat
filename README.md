# astrocat
Testbed PostgreSQL extension for experimenting on indexing astronomical catalogues.

It is heavily based on pgSphere, whose up-to-date code is available at https://github.com/akorotkov/pgsphere

Currently supported are simple radial (a.k.a. *find all points inside given circle*) and kNN ordering (*find given number of closest points to given location*) queries. 

Plans are to implement some sort of spatial clustering optimized for rapid grouping of object measurements in large-scale time-domain surveys, when you have a lot of groups of points with close but not exactly same positions, representing different observations of the same objects at different moments of time. 

### Usage ###

```sql
CREATE EXTENSION astrocat;

-- coordinates are in degreed
CREATE TABLE test (ra FLOAT, dec FLOAT);

-- functional index, does not need any dedicated column with specific type
CREATE INDEX ON test USING gist (skypoint(ra, dec));

-- cone search or radial query
SELECT count(*) FROM test WHERE skypoint(ra,dec) @ skycircle(11.2, 67.0, 2.0);

-- k-nearest-neighbour (knn) ordering
SELECT *, skypoint(ra,dec) <-> skypoint(11.2, 67.0) AS dist FROM test ORDER BY skypoint(ra,dec) <-> skypoint(11.2, 67.0) limit 10;
```
