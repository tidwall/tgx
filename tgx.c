// https://github.com/tidwall/tgx
//
// Copyright 2024 Joshua J Baker. All rights reserved.
// Use of this source code is governed by an MIT-style
// license that can be found in the LICENSE file.
//
// tgx: extensions for the TG library

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tgx.h"

void *tg_malloc(size_t nbytes);
void tg_free(void *ptr);

static GEOSGeometry *tg_point_to_geos(GEOSContextHandle_t handle,
    struct tg_point point)
{
    return GEOSGeom_createPointFromXY_r(handle, point.x, point.y);
}

static GEOSGeometry *tg_ring_to_geos(GEOSContextHandle_t handle,
    const struct tg_ring *ring, bool isline)
{
    unsigned int npoints = (unsigned int)tg_ring_num_points(ring);
    const struct tg_point *points = tg_ring_points(ring);
    GEOSCoordSequence *seq = GEOSCoordSeq_create_r(handle, npoints, 2);
    if (!seq) {
        return 0;
    }
    for (unsigned int i = 0; i < npoints; i++) {
        if (!GEOSCoordSeq_setXY_r(handle, seq, i, points[i].x, points[i].y)) {
            GEOSCoordSeq_destroy_r(handle, seq);
            return 0;
        }
    }
    GEOSGeometry *geom = 0;
    if (isline) {
        geom = GEOSGeom_createLineString_r(handle, seq);
    } else {
        geom = GEOSGeom_createLinearRing_r(handle, seq);
    }
    if (!geom) {
        GEOSCoordSeq_destroy_r(handle, seq);
        return 0;
    }
    return geom;
}

static GEOSGeometry *tg_line_to_geos(GEOSContextHandle_t handle,
    const struct tg_line *line)
{
    return tg_ring_to_geos(handle, (struct tg_ring*)line, true);
}

static GEOSGeometry *tg_poly_to_geos(GEOSContextHandle_t handle,
    const struct tg_poly *poly)
{
    const struct tg_ring *ring = tg_poly_exterior(poly);
    if (!ring) {
        return 0;
    }
    GEOSGeometry *shell = tg_ring_to_geos(handle, ring, false);
    if (!shell) {
        return 0;
    }
    unsigned int nholes = (unsigned int)tg_poly_num_holes(poly);
    GEOSGeometry **holes = 0;
    if (nholes > 0) {
        holes = tg_malloc(nholes * sizeof(GEOSGeometry*));
        if (!holes) {
            GEOSGeom_destroy_r(handle, shell);
            return 0;
        }
        bool ok = true;
        unsigned int i = 0;
        for (; i < nholes; i++) {
            const struct tg_ring *hole = tg_poly_hole_at(poly, i);
            if (!hole) {
                ok = false;
                break;
            }
            holes[i] = tg_ring_to_geos(handle, hole, false);
            if (!holes[i]) {
                ok = false;
                break;
            }
        }
        if (!ok) {
            for (unsigned int j = 0; j < i; j++) {
                GEOSGeom_destroy_r(handle, holes[j]);
            }
            GEOSGeom_destroy_r(handle, shell);
            tg_free(holes);
            return 0;
        }
    }
    GEOSGeometry *geom = GEOSGeom_createPolygon_r(handle, shell, holes, nholes);
    if (!geom) {
        for (unsigned int i = 0; i < nholes; i++) {
            GEOSGeom_destroy_r(handle, holes[i]);
        }
        GEOSGeom_destroy_r(handle, shell);
    }
    if (holes) {
        tg_free(holes);
    }
    return geom;
}

/// Convert a GEOS geometry to a TG geometry
GEOSGeometry *tg_geom_to_geos(GEOSContextHandle_t handle,
    const struct tg_geom *geom)
{
    if (tg_geom_error(geom)) {
        return 0;
    }
    // Attempt to create a simple GEOS geometry without WKB convertion.
    if (!tg_geom_is_empty(geom) && tg_geom_dims(geom) == 2) {
        switch (tg_geom_typeof(geom)) {
        case TG_POINT:
            return tg_point_to_geos(handle, tg_geom_point(geom));
        case TG_LINESTRING:
            return tg_line_to_geos(handle, tg_geom_line(geom));
        case TG_POLYGON:
            return tg_poly_to_geos(handle, tg_geom_poly(geom));
        default:
            break;
        }
    }
    // Everything else needs WKB conversion
    bool must_free = false;
    uint8_t smallfry[4096];
    uint8_t *wkb = smallfry;
    size_t size = tg_geom_wkb(geom, smallfry, sizeof(smallfry));
    if (size > sizeof(smallfry)) {
        uint8_t *wkb2 = tg_malloc(size);
        if (!wkb2) {
            return 0;
        }
        tg_geom_wkb(geom, wkb2, size);
        wkb = wkb2;
        must_free = true;
    }
    GEOSGeometry *ggeom = GEOSGeomFromWKB_buf_r(handle, wkb, size);
    if (must_free) {
        tg_free(wkb);
    }
    return ggeom;
}

static struct tg_geom *tg_point_from_geos(GEOSContextHandle_t handle,
    GEOSGeometry *geom)
{
    double x, y;
    if (GEOSGeomGetX_r(handle, geom, &x) && GEOSGeomGetX_r(handle, geom, &y)) {
        return tg_geom_new_point((struct tg_point){ x, y });
    }
    return 0;
}

static struct tg_ring *tg_ring_from_geos(GEOSContextHandle_t handle,
    const GEOSGeometry *geom, bool isline)
{
    if (!geom) {
        return 0;
    }
    const GEOSCoordSequence *seq = GEOSGeom_getCoordSeq_r(handle, geom);
    if (!seq) {
        return 0;
    }
    unsigned int npoints;
    if (!GEOSCoordSeq_getSize_r(handle, seq, &npoints)) {
        return 0;
    }
    struct tg_point *points = tg_malloc(npoints*sizeof(struct tg_point));
    if (!points) {
        return 0;
    }
    if (!GEOSCoordSeq_copyToBuffer_r(handle, seq, (double*)points, 0, 0)) {
        tg_free(points);
        return 0;
    }
    struct tg_ring *ring = 0;
    if (isline) {
        ring = (struct tg_ring*)tg_line_new(points, (int)npoints);
    } else {
        ring = (struct tg_ring*)tg_ring_new(points, (int)npoints);
    }
    tg_free(points);
    return ring;
}

static struct tg_geom *tg_line_from_geos(GEOSContextHandle_t handle,
    GEOSGeometry *geom)
{
    return (struct tg_geom*)tg_ring_from_geos(handle, geom, true);
}

static struct tg_geom *tg_poly_from_geos(GEOSContextHandle_t handle,
    const GEOSGeometry *geom)
{
    struct tg_ring *exterior = tg_ring_from_geos(handle, 
        GEOSGetExteriorRing_r(handle, geom), false);
    if (!exterior) {
        return 0;
    }
    int nholes = GEOSGetNumInteriorRings_r(handle, geom);
    if (nholes < 0) {
        tg_ring_free(exterior);
        return 0;
    }
    struct tg_ring **holes = tg_malloc(nholes*sizeof(struct tg_ring*));
    if (!holes) {
        tg_ring_free(exterior);
        return 0;
    }
    int i = 0;
    for (; i < nholes; i++) {
        holes[i] = tg_ring_from_geos(handle,
            GEOSGetInteriorRingN_r(handle, geom, i), false);
        if (!holes[i]) {
            break;
        }
    }
    struct tg_poly *poly = 0;
    if (i == nholes) {
        poly = tg_poly_new(exterior, (void*)holes, nholes);
    }
    for (int j = 0; j < i; j++) {
        tg_ring_free(holes[j]);
    }
    tg_free(holes);
    tg_ring_free(exterior);
    return (struct tg_geom*)poly;
}

/// Convert a TG geometry to a GEOS geometry
struct tg_geom *tg_geom_from_geos(GEOSContextHandle_t handle,
    GEOSGeometry *geom)
{
    if (!geom) {
        return 0;
    }
    if (!GEOSisEmpty_r(handle, geom) && 
        GEOSGeom_getDimensions_r(handle, geom) == 2)
    {
        switch (GEOSGeomTypeId_r(handle, geom)) {
        case GEOS_POINT:
            return tg_point_from_geos(handle, geom);
        case GEOS_LINESTRING:
            return tg_line_from_geos(handle, geom);
        case GEOS_POLYGON:
            return tg_poly_from_geos(handle, geom);
        default:
            break;
        }
    }
    // Everything else needs WKB conversion
    size_t size;
    unsigned char *wkb = GEOSGeomToWKB_buf_r(handle, geom, &size);
    if (!wkb) {
        return 0;
    }
    struct tg_geom *tgeom = tg_parse_wkb((uint8_t*)wkb, size);
    GEOSFree_r(handle, wkb);
    return tgeom;
}

#define EARTH 6378137                // WGS84 equatorial radius in meters
#define FLATTENING (1/298.257223563) // WGS84 flattening

#define RADIANS (M_PI / 180)
#define DEGREES (180 / M_PI)

static double sphere_dist(double lat1, double lon1, double lat2, double lon2) {
    // haversine
    double rlat1 = lat1 * RADIANS;
    double rlat2 = lat2 * RADIANS;
    double rlat = (lat2-lat1) * RADIANS;
    double rlon = (lon2-lon1) * RADIANS;
    double a = sin(rlat/2) * sin(rlat/2) + 
        cos(rlat1) * cos(rlat2) * sin(rlon/2) * sin(rlon/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return EARTH * c;
}

static double sphere_initial_azi(double lat1, double lon1, double lat2,
    double lon2)
{
    double rlat1 = lat1 * RADIANS;
    double rlat2 = lat2 * RADIANS;
    double rlon = (lon2 - lon1) * RADIANS;
    double y = sin(rlon) * cos(rlat2);
    double x = cos(rlat1) * sin(rlat2) - sin(rlat1) * cos(rlat2) * cos(rlon);
    double a = atan2(y, x);
    double bearing = a * DEGREES;
    return fmod(bearing+360, 360);
}

static void sphere_dest(double lat1, double lon1, double azi1, double s12,
    double *lat2, double *lon2)
{
    double dr = s12 / EARTH;
    double ar = azi1 * RADIANS;
    double rlat1 = lat1 * RADIANS;
    double rlon1 = lon1 * RADIANS;
    double phi = asin(sin(rlat1) * cos(dr) + cos(rlat1) * sin(dr) * cos(ar));
    double lambda = rlon1 + atan2(sin(ar) * sin(dr) * cos(rlat1), cos(dr) - 
        sin(rlat1) * sin(phi));
    lambda = fmod(lambda+3*M_PI, 2*M_PI) - M_PI; // normalise to -180..+180°
    *lat2 = phi * DEGREES;
    *lon2 = lambda * DEGREES;
}

static double geo_dist(struct tg_point a, struct tg_point b) {
    return sphere_dist(a.y, a.x, b.y, b.x);
}

static double geo_bear(struct tg_point a, struct tg_point b) {
    return sphere_initial_azi(a.y, a.x, b.y, b.x);
}

static struct tg_point geo_dest(struct tg_point origin, double bearing,
    double meters)
{
    double y, x;
    sphere_dest(origin.y, origin.x, bearing, meters, &y, &x);
    return (struct tg_point){ x, y };
}

static bool feq(double x, double y) {
    return !((x < y) | (x > y));
}

static bool pteq(struct tg_point a, struct tg_point b) {
    return feq(a.x, b.x) && feq(a.y, b.y);
}

static struct tg_point point_to_meters_grid(struct tg_point point,
    struct tg_point origin)
{
    if (pteq(point, origin)) {
        return (struct tg_point) { 0, 0 };
    }
    double dist = geo_dist(origin, point);
    double bear = geo_bear(origin, point);
    double x = dist * sin(bear * RADIANS);
    double y = dist * cos(bear * RADIANS);
    return (struct tg_point) { x, y };
}

static struct tg_point point_from_meters_grid(struct tg_point point,
    struct tg_point origin)
{
    struct tg_point result = { 0, 0 };    
    double x = point.x;
    double y = point.y;
    double bear = atan2(x, y)*DEGREES;
    // bear = fmod(bear+360, 360); // TODO << is this needed?
    double dist = sqrt(x * x + y * y);
    result = geo_dest(origin, bear, dist);
    return result;
}

static struct tg_geom *tg_point_to_meters_grid(struct tg_point point, 
    struct tg_point origin)
{
    return tg_geom_new_point(point_to_meters_grid(point, origin));
}

static struct tg_ring *ring_to_meters_grid(const struct tg_ring *ring, 
    struct tg_point origin, bool isline)
{
    int npoints = tg_ring_num_points(ring);
    const struct tg_point *points = tg_ring_points(ring);
    struct tg_point *points2 = tg_malloc(npoints*sizeof(struct tg_point));
    if (!points2) {
        return 0;
    }
    for (int i = 0; i < npoints; i++) {
        points2[i] = point_to_meters_grid(points[i], origin);
    }
    struct tg_ring *ring2 = 0;
    if (isline) {
        ring2 = (struct tg_ring*)tg_line_new(points2, npoints);
    } else {
        ring2 = (struct tg_ring*)tg_ring_new(points2, npoints);
    }
    tg_free(points2);
    return ring2;
}

static struct tg_line *line_to_meters_grid(const struct tg_line *line,
    struct tg_point origin)
{
    return (struct tg_line*)ring_to_meters_grid((struct tg_ring*)line, origin, 
        true);
}

static struct tg_geom *tg_line_to_meters_grid(const struct tg_line *line,
    struct tg_point origin)
{
    return (struct tg_geom*)line_to_meters_grid(line, origin);
}

static struct tg_poly *poly_to_meters_grid(const struct tg_poly *poly,
    struct tg_point origin)
{
    struct tg_ring *exterior = ring_to_meters_grid(tg_poly_exterior(poly), 
        origin, false);
    if (!exterior) {
        return 0;
    }
    int nholes = tg_poly_num_holes(poly);
    struct tg_ring **holes = tg_malloc(nholes*sizeof(struct tg_ring*));
    if (!holes) {
        tg_ring_free(exterior);
        return 0;
    }
    int i = 0;
    for (; i < nholes; i++) {
        holes[i] = ring_to_meters_grid(tg_poly_hole_at(poly, i), origin, false);
        if (!holes[i]) {
            break;
        }
    }
    struct tg_poly *poly2 = 0;
    if (i == nholes) {
        poly2 = tg_poly_new(exterior, (void*)holes, nholes); 
    }
    for (int j = 0; j < i; j++) {
        tg_ring_free(holes[j]);
    }
    tg_free(holes);
    tg_ring_free(exterior);
    return poly2;
}

static struct tg_geom *tg_poly_to_meters_grid(const struct tg_poly *poly,
    struct tg_point origin)
{
    return (struct tg_geom*)poly_to_meters_grid(poly, origin);
}

static struct tg_geom *tg_multipoint_to_meters_grid(const struct tg_geom *geom,
    struct tg_point origin)
{
    int npoints = tg_geom_num_points(geom);
    struct tg_point *points = tg_malloc(npoints*sizeof(struct tg_point));
    if (!points) {
        return 0;
    }
    for (int i = 0; i < npoints; i++) {
        points[i] = point_to_meters_grid(tg_geom_point_at(geom, i), origin);
    }
    struct tg_geom *geom2 = tg_geom_new_multipoint(points, npoints);
    tg_free(points);
    return geom2;
}

static struct tg_geom *tg_multilinestring_to_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int nlines = tg_geom_num_lines(geom);
    struct tg_line **lines = tg_malloc(nlines*sizeof(struct tg_line*));
    if (!lines) {
        return 0;
    }
    int i = 0;
    for (; i < nlines; i++) {
        lines[i] = line_to_meters_grid(tg_geom_line_at(geom, i), origin);
        if (!lines[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == nlines) {
        geom2 = tg_geom_new_multilinestring((void*)lines, nlines);
    }
    for (int j = 0; j < i; j++) {
        tg_line_free(lines[j]);
    }
    tg_free(lines);
    return geom2;
}

static struct tg_geom *tg_multipolygon_to_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int npolys = tg_geom_num_polys(geom);
    struct tg_poly **polys = tg_malloc(npolys*sizeof(struct tg_poly*));
    if (!polys) {
        return 0;
    }
    int i = 0;
    for (; i < npolys; i++) {
        polys[i] = poly_to_meters_grid(tg_geom_poly_at(geom, i), origin);
        if (!polys[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == npolys) {
        geom2 = tg_geom_new_multipolygon((void*)polys, npolys);
    }
    for (int j = 0; j < i; j++) {
        tg_poly_free(polys[j]);
    }
    tg_free(polys);
    return geom2;
}

static struct tg_geom *tg_geometrycollection_to_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int ngeoms = tg_geom_num_geometries(geom);
    struct tg_geom **geoms = tg_malloc(ngeoms*sizeof(struct tg_geom*));
    if (!geoms) {
        return 0;
    }
    int i = 0;
    for (; i < ngeoms; i++) {
        geoms[i] = tg_geom_to_meters_grid(tg_geom_geometry_at(geom, i),
            origin);
        if (!geoms[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == ngeoms) {
        geom2 = tg_geom_new_geometrycollection((void*)geoms, ngeoms);
    }
    for (int j = 0; j < i; j++) {
        tg_geom_free(geoms[j]);
    }
    tg_free(geoms);
    return geom2;
}

/// Convert a TG geometry from WGS84 coordinates to meters grid
/// @param geom the geometry to convert
/// @param origin a WGS84 origin point for the input geometry
struct tg_geom *tg_geom_to_meters_grid(const struct tg_geom *geom,
    struct tg_point origin)
{
    if (!tg_geom_error(geom)) {
        switch (tg_geom_typeof(geom)) {
        case TG_POINT:
            return tg_point_to_meters_grid(tg_geom_point(geom), origin);
        case TG_LINESTRING:
            return tg_line_to_meters_grid(tg_geom_line(geom), origin);
        case TG_POLYGON:
            return tg_poly_to_meters_grid(tg_geom_poly(geom), origin);
        case TG_MULTIPOINT:
            return tg_multipoint_to_meters_grid(geom, origin);
        case TG_MULTILINESTRING:
            return tg_multilinestring_to_meters_grid(geom, origin);
        case TG_MULTIPOLYGON:
            return tg_multipolygon_to_meters_grid(geom, origin);
        case TG_GEOMETRYCOLLECTION:
            return tg_geometrycollection_to_meters_grid(geom, origin);
        }
    }
    return 0;
}








////////////////////////////////


static struct tg_geom *tg_point_from_meters_grid(struct tg_point point, 
    struct tg_point origin)
{
    return tg_geom_new_point(point_from_meters_grid(point, origin));
}

static struct tg_ring *ring_from_meters_grid(const struct tg_ring *ring, 
    struct tg_point origin, bool isline)
{
    int npoints = tg_ring_num_points(ring);
    const struct tg_point *points = tg_ring_points(ring);
    struct tg_point *points2 = tg_malloc(npoints*sizeof(struct tg_point));
    if (!points2) {
        return 0;
    }
    for (int i = 0; i < npoints; i++) {
        points2[i] = point_from_meters_grid(points[i], origin);
    }
    struct tg_ring *ring2 = 0;
    if (isline) {
        ring2 = (struct tg_ring*)tg_line_new(points2, npoints);
    } else {
        ring2 = (struct tg_ring*)tg_ring_new(points2, npoints);
    }
    tg_free(points2);
    return ring2;
}

static struct tg_line *line_from_meters_grid(const struct tg_line *line,
    struct tg_point origin)
{
    return (struct tg_line*)ring_from_meters_grid((struct tg_ring*)line, origin, 
        true);
}

static struct tg_geom *tg_line_from_meters_grid(const struct tg_line *line,
    struct tg_point origin)
{
    return (struct tg_geom*)line_from_meters_grid(line, origin);
}

static struct tg_poly *poly_from_meters_grid(const struct tg_poly *poly,
    struct tg_point origin)
{
    struct tg_ring *exterior = ring_from_meters_grid(tg_poly_exterior(poly), 
        origin, false);
    if (!exterior) {
        return 0;
    }
    int nholes = tg_poly_num_holes(poly);
    struct tg_ring **holes = tg_malloc(nholes*sizeof(struct tg_ring*));
    if (!holes) {
        tg_ring_free(exterior);
        return 0;
    }
    int i = 0;
    for (; i < nholes; i++) {
        holes[i] = ring_from_meters_grid(tg_poly_hole_at(poly, i), origin,
            false);
        if (!holes[i]) {
            break;
        }
    }
    struct tg_poly *poly2 = 0;
    if (i == nholes) {
        poly2 = tg_poly_new(exterior, (void*)holes, nholes); 
    }
    for (int j = 0; j < i; j++) {
        tg_ring_free(holes[j]);
    }
    tg_free(holes);
    tg_ring_free(exterior);
    return poly2;
}

static struct tg_geom *tg_poly_from_meters_grid(const struct tg_poly *poly,
    struct tg_point origin)
{
    return (struct tg_geom*)poly_from_meters_grid(poly, origin);
}

static struct tg_geom *tg_multipoint_from_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int npoints = tg_geom_num_points(geom);
    struct tg_point *points = tg_malloc(npoints*sizeof(struct tg_point));
    if (!points) {
        return 0;
    }
    for (int i = 0; i < npoints; i++) {
        points[i] = point_from_meters_grid(tg_geom_point_at(geom, i), origin);
    }
    struct tg_geom *geom2 = tg_geom_new_multipoint(points, npoints);
    tg_free(points);
    return geom2;
}

static struct tg_geom *tg_multilinestring_from_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int nlines = tg_geom_num_lines(geom);
    struct tg_line **lines = tg_malloc(nlines*sizeof(struct tg_line*));
    if (!lines) {
        return 0;
    }
    int i = 0;
    for (; i < nlines; i++) {
        lines[i] = line_from_meters_grid(tg_geom_line_at(geom, i), origin);
        if (!lines[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == nlines) {
        geom2 = tg_geom_new_multilinestring((void*)lines, nlines);
    }
    for (int j = 0; j < i; j++) {
        tg_line_free(lines[j]);
    }
    tg_free(lines);
    return geom2;
}

static struct tg_geom *tg_multipolygon_from_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int npolys = tg_geom_num_polys(geom);
    struct tg_poly **polys = tg_malloc(npolys*sizeof(struct tg_poly*));
    if (!polys) {
        return 0;
    }
    int i = 0;
    for (; i < npolys; i++) {
        polys[i] = poly_from_meters_grid(tg_geom_poly_at(geom, i), origin);
        if (!polys[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == npolys) {
        geom2 = tg_geom_new_multipolygon((void*)polys, npolys);
    }
    for (int j = 0; j < i; j++) {
        tg_poly_free(polys[j]);
    }
    tg_free(polys);
    return geom2;
}

static struct tg_geom *tg_geometrycollection_from_meters_grid(
    const struct tg_geom *geom, struct tg_point origin)
{
    int ngeoms = tg_geom_num_geometries(geom);
    struct tg_geom **geoms = tg_malloc(ngeoms*sizeof(struct tg_geom*));
    if (!geoms) {
        return 0;
    }
    int i = 0;
    for (; i < ngeoms; i++) {
        geoms[i] = tg_geom_from_meters_grid(tg_geom_geometry_at(geom, i),
            origin);
        if (!geoms[i]) {
            break;
        }
    }
    struct tg_geom *geom2 = 0;
    if (i == ngeoms) {
        geom2 = tg_geom_new_geometrycollection((void*)geoms, ngeoms);
    }
    for (int j = 0; j < i; j++) {
        tg_geom_free(geoms[j]);
    }
    tg_free(geoms);
    return geom2;
}

/// Convert a TG geometry from meters grid to WGS84 coordinates
/// @param geom the geometry to convert
/// @param origin a WGS84 origin point for the output geometry
struct tg_geom *tg_geom_from_meters_grid(const struct tg_geom *geom,
    struct tg_point origin)
{
    if (!tg_geom_error(geom)) {
        switch (tg_geom_typeof(geom)) {
        case TG_POINT:
            return tg_point_from_meters_grid(tg_geom_point(geom), origin);
        case TG_LINESTRING:
            return tg_line_from_meters_grid(tg_geom_line(geom), origin);
        case TG_POLYGON:
            return tg_poly_from_meters_grid(tg_geom_poly(geom), origin);
        case TG_MULTIPOINT:
            return tg_multipoint_from_meters_grid(geom, origin);
        case TG_MULTILINESTRING:
            return tg_multilinestring_from_meters_grid(geom, origin);
        case TG_MULTIPOLYGON:
            return tg_multipolygon_from_meters_grid(geom, origin);
        case TG_GEOMETRYCOLLECTION:
            return tg_geometrycollection_from_meters_grid(geom, origin);
        }
    }
    return 0;
}
