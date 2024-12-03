# tgx

Extention for the [TG](https://github.com/tidwall/tg) geometry library.

Includes a bridge for converting geometries to/from the [GEOS](https://github.com/libgeos/geos) library and simple a utility for projecting geometries from lat/lon (spherical)
to cartesian coordinates using a meters grid.

Required TG and GEOS to use.

## Example

Here's an example of how to take a GeoJSON geometry of the city of Tempe, 
Arizona and buffer its boundary by 1000 meters.


```c
// Parse the geojson and convert from lat/lon to meters.
struct tg_geom *tg_geom = tg_parse_geojson(tempe);
struct tg_point origin = tg_geom_point(tg_geom);
struct tg_geom *tg_meters = tg_geom_to_meters_grid(tg_geom, origin);

// Convert TG geometry to GEOS geometry
GEOSContextHandle_t handle = GEOS_init_r();
GEOSGeometry *geos_geom = tg_geom_to_geos(handle, tg_meters);

// Buffer the geometry using the GEOS buffer operation. 
// The width is 1000 meters.
GEOSGeometry *geos_buffered = GEOSBuffer_r(handle, geos_geom, 1000, 20);

// Convert back to TG from GEOS
struct tg_geom *tg_buffered = tg_geom_from_geos(handle, geos_buffered);

// Convert back to a lat/lon geometry.
// This is the final result.
struct tg_geom *tg_final = tg_geom_from_meters_grid(tg_buffered, origin);

// Finally, cleanup resources
GEOSGeom_destroy_r(handle, geos_buffered);
GEOSGeom_destroy_r(handle, geos_geom);
GEOS_finish_r(handle);

tg_geom_free(tg_buffered);
tg_geom_free(tg_meters);
tg_geom_free(tg_geom);
tg_geom_free(tg_final);
```

<img width="635" alt="image" src="https://github.com/user-attachments/assets/e522b252-4106-4c9e-881f-b236accf035c">


