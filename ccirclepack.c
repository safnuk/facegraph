// ccirclepack.c

#include <stdlib.h>
#include <stdio.h>
#include "cmesh.h"
#include "ccirclepack.h"

float calc_circlepack_metric(mesh *m)
{
        float *ambient_lengths;
        ambient_lengths = (float *)malloc(m->ranks[1] *sizeof(float));
        if (ambient_lengths == NULL) {
                printf("Memory allocation error");
                exit(1);
        }
        calc_ambient_edge_lengths(m, ambient_lengths);
        for (int i=0; i < m->ranks[1]; i++) {
                printf("%f, ", ambient_lengths[i]);
        }
        free(ambient_lengths)
}

void calc_ambient_edge_lengths(mesh *m, float *edge_lengths)
{
        int i;
        point *v1, *v2;
        for (i=0; i < m->ranks[1]; i++)
        {
                v1 = get_coordinate(m, (m->edges[i]).vertices[0]);
                v2 = get_coordinate(m, (m->edges[i]).vertices[1]);
                edge_lengths[i] = calc_distance(v1, v2);
        }
}

void calc_initial_radii(mesh *m, float *edge_lengths)
{
}
