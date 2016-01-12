MODULES = astrocat

EXTENSION = astrocat
DATA = astrocat--1.0.sql
PGFILEDESC = "astrocat - astronomical catalogue extension for postgresql"

REGRESS = astrocat

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
