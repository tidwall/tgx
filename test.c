// cc -Ideps -Ideps/libgeos/build/install/include *.c deps/*.c deps/libgeos/build/install/lib/libgeos_c.a deps/libgeos/build/install/lib/libgeos.a -lm  -lstdc++ && ./a.out

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
#include <string.h>
#include <sys/types.h>
#include <errno.h>
#include "tgx.h"

static int tmallocs = 0;

void *xmalloc(size_t size) {
    tmallocs++;
    return malloc(size);
}

void *xrealloc(void *ptr, size_t size) {
    if (!ptr) {
        return xmalloc(size);
    }
    return realloc(ptr, size);
}

void xfree(void *ptr) {
    if (ptr) {
        tmallocs--;
        free(ptr);
    }
}

static void test_geojson_file(const char *name, const char *gstr) {
    (void)name;
    GEOSContextHandle_t handle = GEOS_init_r();
    struct tg_geom *tgeom = tg_parse_geojson(gstr);
    assert(!tg_geom_error(tgeom));
    GEOSGeometry *ggeom = tg_geom_to_geos(handle, tgeom);
    tg_geom_free(tgeom);
    assert(ggeom);

    GEOSGeometry *ggeom2 = GEOSBuffer_r(handle, ggeom, 10, 10);
    GEOSGeom_destroy_r(handle, ggeom);
    assert(ggeom2);
    struct tg_geom *tgeom2 = tg_geom_from_geos(handle, ggeom2);
    GEOSGeom_destroy_r(handle, ggeom2);
    GEOS_finish_r(handle);
    assert(!tg_geom_error(tgeom2));

    char output[4096];
    tg_geom_geojson(tgeom2, output, sizeof(output));


    struct tg_point origin = tg_geom_point(tgeom2);
    struct tg_geom *tgeom3 = tg_geom_to_meters_grid(tgeom2, origin);
    tg_geom_free(tgeom2);
    assert(tgeom3);

    struct tg_geom *tgeom4 = tg_geom_from_meters_grid(tgeom3, origin);
    tg_geom_free(tgeom3);
    assert(tgeom4);


    tg_geom_free(tgeom4);

}

static char *read_testfile(const char *name) {
    char path[1024];
    snprintf(path, sizeof(path), "testfiles/%s", name);
    FILE *f = fopen(path, "rb+");
    if (!f) {
        fprintf(stderr, "%s: %s\n", strerror(errno), name);
        exit(1);
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    rewind(f);
    char *data = malloc(size+1);
    assert(data);
    fread(data, 1, size, f);
    data[size] = 0;
    fclose(f);
    return data;
}

static void test_geojson_dir(void) {
    struct dirent *entry;
    DIR *dir = opendir("testfiles");
    if (!dir) {
        fprintf(stderr, "%s\n", strerror(errno));
        exit(1);
    }
    while ((entry = readdir(dir))) {
        if (entry->d_name[0] == '.') {
            continue;
        }
        char *data = read_testfile(entry->d_name);
        test_geojson_file(entry->d_name, data);
        free(data);
    }
    closedir(dir);

    assert(tmallocs == 0);
}


void print_geojson(struct tg_geom *geom) {
    size_t n = tg_geom_geojson(geom, 0, 0);
    char *output = malloc(n+1);
    assert(output);
    tg_geom_geojson(geom, output, n+1);
    printf("%s\n", output);
    free(output);
}

int main(int argc, char *argv[]) {
    (void)argc, (void)argv;
    tg_env_set_allocator(xmalloc, xrealloc, xfree);
    test_geojson_dir();





    // char *gstr = read_testfile("tempe.geojson");
    // // char *gstr = 
    // // "{\"type\":\"Polygon\",\"coordinates\":[["
    // // "[-112.5954,32.7146],[-112.5954,32.7143],[-112.5949,32.7143],"
    // // "[-112.5949,32.7146],[-112.5954,32.7146]"
    // // "]]}";
    // // "{\"type\":\"LineString\",\"coordinates\":["
    // // "[-112.5954,32.7146],[-112.5949,32.7143]"
    // // "]}";
    // // "{\"type\":\"Point\",\"coordinates\":[-112.5954,32.7146]}";
    

    // double meters = 10;

    // struct tg_geom *geom = tg_parse(gstr, strlen(gstr));
    // assert(!tg_geom_error(geom));
    // // free(gstr);

    // // {
    // //     struct tg_point center = tg_geom_point(geom);
    // //     struct tg_geom *gm = tg_geom_to_meters_grid(geom, center);
    // //     // struct tg_geom *gr = tg_geom_buffer(gm, 10, 20);
    // //     struct tg_geom *gg = tg_geom_from_meters_grid(gr, center);


        
    // //     // 
    // //     print_geojson(gg);
    // //     exit(0);
    // // }

    // struct tg_point origin = tg_geom_point(geom);
    // // struct tg_point origin = { -112.59515,32.71445 };

    // struct tg_geom *geom2 = tg_geom_to_meters_grid(geom, origin);
    // assert(!tg_geom_error(geom2));

    // GEOSContextHandle_t handle = GEOS_init_r();
    // GEOSGeometry *ggeom = tg_geom_to_geos(handle, geom2);
    // assert(ggeom);
    // GEOSGeometry *ggeom2 = GEOSBuffer_r(handle, ggeom, 1000, 20);
    // GEOSGeom_destroy_r(handle, ggeom);
    // assert(ggeom2);
    // struct tg_geom *geom3 = tg_geom_from_geos(handle, ggeom2);
    // GEOSGeom_destroy_r(handle, ggeom2);
    // GEOS_finish_r(handle);
    // assert(!tg_geom_error(geom3));
    


    // struct tg_geom *geom4 = tg_geom_from_meters_grid(geom3, origin);
    // assert(!tg_geom_error(geom4));
    // // printf("{\"type\":\"GeometryCollection\",\"geometries\":[\n");
    // // print_geojson(geom);
    // // printf(",");
    // // print_geojson(tg_geom_new_point(origin));
    // // printf(",");
    // // print_geojson(geom4);
    // // printf("]}\n");
    


    // {
    // // Parse the geojson and convert from lat/lon to meters.
    // struct tg_geom *tg_geom = tg_parse_geojson(gstr);
    // struct tg_point origin = tg_geom_point(tg_geom);
    // struct tg_geom *tg_meters = tg_geom_to_meters_grid(tg_geom, origin);

    // // Convert TG geometry to GEOS geometry
    // GEOSContextHandle_t handle = GEOS_init_r();
    // GEOSGeometry *geos_geom = tg_geom_to_geos(handle, tg_meters);

    // // Buffer the geometry using the GEOS buffer operation. 
    // // The width is 1000 meters.
    // GEOSGeometry *geos_buffered = GEOSBuffer_r(handle, geos_geom, 1000, 20);

    // // Convert back to TG from GEOS
    // struct tg_geom *tg_buffered = tg_geom_from_geos(handle, geos_buffered);

    // // Convert back to lat/lon.
    // // This is the final result.
    // struct tg_geom *tg_final = tg_geom_from_meters_grid(tg_buffered, origin);



    // // Finally, cleanup resources
    // GEOSGeom_destroy_r(handle, geos_buffered);
    // GEOSGeom_destroy_r(handle, geos_geom);
    // GEOS_finish_r(handle);
    
    // tg_geom_free(tg_buffered);
    // tg_geom_free(tg_meters);
    // tg_geom_free(tg_geom);
    // tg_geom_free(tg_final);



    // }


    // // // tg_parse_geojson
    // // GEOSContextHandle_t handle = GEOS_init_r();
    // // char *gstr = read_testfile("point.geojson");
    // // struct tg_geom *tgeom = tg_parse(gstr, strlen(gstr));
    // // free(gstr);
    // // assert(!tg_geom_error(tgeom));
    // // GEOSGeometry *ggeom = tg_to_geos(handle, tgeom);
    // // tg_geom_free(tgeom);
    // // assert(ggeom);

    // // GEOSGeometry *ggeom2 = GEOSBuffer_r(handle, ggeom, 10, 10);
    // // GEOSGeom_destroy_r(handle, ggeom);
    // // assert(ggeom2);
    // // struct tg_geom *tgeom2 = tg_from_geos(handle, ggeom2);
    // // GEOSGeom_destroy_r(handle, ggeom2);
    // // GEOS_finish_r(handle);
    // // assert(!tg_geom_error(tgeom2));

    // // char output[4096];
    // // tg_geom_geojson(tgeom2, output, sizeof(output));
    // // tg_geom_free(tgeom2);
    // // printf("%s\n", output);

    return 0;
}