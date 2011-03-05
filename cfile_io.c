/* cfile_io.c
 *
*/

#include <stdio.h>

typedef struct _point {
        float x, y, z;
        struct _point *next;
} point;

typedef struct _face {
        int v1, v2, v3;
        struct _face *next;
} face;

typedef struct {
        int number_of_points;
        int number_of_triangles;
        float *coordinates;
        int *faces;
} filedata; 

int read(char *filename, filedata *data)
{
        FILE *fp;
        int point_count = 0;
        int face_count = 0;

        fp = fopen(filename, "r");
        if (fp == NULL) {
                fprintf(stderr, "Can't open file!\n");
                exit(1);
        }

}
