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

// parse a wkt using TG, convert to geos and back to TG, perform wkbs of before
// after geos<->tg conversion
static void proc_conv(const char *wkt) {
    // printf("%s\n", wkt);
    uint8_t wkb1[512];
    uint8_t wkb2[512];
    struct tg_geom *g1 = tg_parse(wkt, strlen(wkt));
    // printf("   A %s\n", tg_geom_error(g1));
    assert(!tg_geom_error(g1));
    size_t n1 = tg_geom_wkb(g1, wkb1, sizeof(wkb1));
    GEOSContextHandle_t handle = GEOS_init_r();
    assert(handle);
    GEOSGeometry *gg = tg_geom_to_geos(handle, g1);
    assert(gg);
    struct tg_geom *g2 = tg_geom_from_geos(handle, gg);
    assert(!tg_geom_error(g2));
    GEOS_finish_r(handle);
    size_t n2 = tg_geom_wkb(g2, wkb2, sizeof(wkb2));
    tg_geom_free(g2);
    tg_geom_free(g1);
    assert(n1 <= sizeof(wkb1));
    assert(n2 <= sizeof(wkb2));
    assert(n1 == n2);
    assert(memcmp(wkb1, wkb2, n1) == 0);
}

static void test_various_conversions(void) {
    proc_conv("POINT(1 2)");
    proc_conv("POINT EMPTY");
    proc_conv("LINESTRING EMPTY");
    proc_conv("LINESTRING(1 2,3 4)");
    proc_conv("POLYGON EMPTY");
    proc_conv("POLYGON((1 2,3 4,5 6,1 2))");
    proc_conv("POLYGON((10 20,30 40,50 60,10 20),(1 2,3 4,5 6,1 2))");
    proc_conv("MULTIPOINT EMPTY");
    proc_conv("MULTIPOINT(1 2,3 4)");
    proc_conv("MULTILINESTRING EMPTY");
    proc_conv("MULTILINESTRING((1 2,3 4),(5 6,7 8))");
    proc_conv("MULTIPOLYGON EMPTY");
    proc_conv("MULTIPOLYGON(((1 2,3 4,5 6,1 2)))");
    proc_conv("MULTIPOLYGON(((10 20,30 40,50 60,10 20),(1 2,3 4,5 6,1 2)),((1 2,3 4,5 6,1 2)))");
    proc_conv("GEOMETRYCOLLECTION EMPTY");
    proc_conv("GEOMETRYCOLLECTION(POINT(1 2))");
    proc_conv("GEOMETRYCOLLECTION(POINT(1 2),MULTIPOINT(1 2,3 4),MULTIPOLYGON(((1 2,3 4,5 6,1 2))))");
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
    test_various_conversions();
    return 0;
}