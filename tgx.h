// https://github.com/tidwall/tgx
//
// Copyright 2024 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.
//
// tgx: extensions for the TG library

#ifndef TGX_H
#define TGX_H

#include <tg.h>
#include <geos_c.h>

GEOSGeometry *tg_geom_to_geos(GEOSContextHandle_t handle, const struct tg_geom *geom);
struct tg_geom *tg_geom_from_geos(GEOSContextHandle_t handle, GEOSGeometry *geom);

struct tg_geom *tg_geom_to_meters_grid(const struct tg_geom *geom, struct tg_point origin);
struct tg_geom *tg_geom_from_meters_grid(const struct tg_geom *geom, struct tg_point origin);

#endif
