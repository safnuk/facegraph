/* cfile_io.c
 *
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>

#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"

#include "cfile_io.h"
#include "geodesic.h"
#include "cmesh.h"



/* Read face data from filename, and store it in data.
 * Assumes that data has been properly initialized.
 */
int read(char *filename, filedata* data) 
{
        char line[100];
        FILE *fp;
        int point_count = 0;
        int face_count = 0;
        _point *point_node = data->point_head;
        face *face_node = data->face_head;

        fp = fopen(filename, "r");
        if (fp == NULL) {
                printf("Can't open file!\n");
                exit(1);
        }
        while (fgets(line, 100, fp)) {
                if (line[0] == 'V') {
                        point_node = parse_new_point(point_node, 
                                        line, point_count);
                        point_count++;
                } else if (line[0] == 'F') {
                        face_node = parse_new_face(face_node, 
                                        line, face_count);
                        face_count++;
                }
        }
        fclose(fp);
        data->number_of_points = point_count;
        data->number_of_triangles = face_count;
        return 0;
}

void save_mesh(char *filename, void *data)
{
        FILE *fp;
        int i, j;
        double r, s;
        int k[3];
        fp = fopen(filename, "w");
        vertex *v;
        edge *e;
        mesh *m = (mesh *) data;
        if (fp == NULL) {
                printf("Can't open file!\n");
                exit(1);
        }
        for (i=0; i < m->ranks[0]; i++) {
                s = m->vertices[i].s;
                r = log((1 + s) / (1 - s));
                fprintf(fp, "Vertex %i %f %f %f %.12f\n", i, m->coordinates[i].x,
                               m->coordinates[i].y, m->coordinates[i].z, r);
        }
        for (i=0; i < m->ranks[2]; i++) {
                for (j=0; j<3; j++) {
                        v = (vertex *)(m->triangles[i].vertices[j]);
                        k[j] = v->index;
                }
                fprintf(fp, "Face %i %i %i %i\n", i, k[0], k[1], k[2]);
        }
        for (i=0; i < m->ranks[1]; i++) {
                e = &(m->edges[i]);
                for (j=0; j<2; j++) {
                        v = (vertex *)(e->vertices[j]);
                        k[j] = v->index;
                }
                if (k[0] > k[1]) {
                        k[2] = k[0];
                        k[0] = k[1];
                        k[1] = k[2];
                }
                fprintf(fp, "Edge %i %i %.12f %.12f\n", k[0], k[1], acos(e->cos_angle),
                      acosh(e->cosh_length) );
        }
        fclose(fp);
}
/* Allocate a new point and add to tail of linked list.
 */
_point *add_point_node(_point *node, double x, double y, double z)
{
        _point *new_node;
        new_node = (_point*)malloc(sizeof(_point));
        if (new_node == NULL) {
                printf("Could not allocate memory for new point.\n");
                exit(1);
        }
        node->next = new_node;
        new_node->next = NULL;
        new_node->x = x;
        new_node->y = y;
        new_node->z = z;
        return new_node;
}

/* Allocate a new face and add to tail of linked list.
 */
face *add_face_node(face *node, int v1, int v2, int v3)
{
        face *new_node;
        new_node = (face*)malloc(sizeof(face));
        if (new_node == NULL) {
                printf("Could not allocate memory for new face.\n");
                exit(1);
        }
        node->next = new_node;
        new_node->next = NULL;
        new_node->v[0] = v1;
        new_node->v[1] = v2;
        new_node->v[2] = v3;
        return new_node;
}

/* Adds a new node to linked list with tail point_node.
 * line should be a string of the form
 *      "Vertex index x y z [radius]"
 * The optional radius value is from a circle packing metric.
 */
_point *parse_new_point(_point* point_node, char *line, int point_count) 
{
        double x, y, z;
        float tx, ty, tz, tr;
        int index;

        sscanf(line, "Vertex %i %f %f %f %f",
                        &index, &tx, &ty, &tz, &tr);
        if (point_count != index) {
                printf("Index mismatch when reading vertex %i", index);
                exit(1);
        }
        x = (double)tx;
        y = (double)ty;
        z = (double)tz;
        return add_point_node(point_node, x, y, z);
}

/* Adds a new node to linked list with tail face_node.
 * line should be a string of the form:
 * "Face index v1 v2 v3"
 * where v1, v2, v3 is a triple of indices to points.
 */
face *parse_new_face(face *face_node, char *line, int face_count)
{
        int v1, v2, v3, index;

        sscanf(line, "Face %i %i %i %i", 
                        &index, &v1, &v2, &v3);
        if (face_count != index) {
                printf("Index mismatch when reading face %i", index);
                exit(1);
        }
        return add_face_node(face_node, v1, v2, v3);
}

void initialize_filedata(filedata *data)
{
        data->point_head = (_point*)malloc(sizeof(_point));
        data->face_head = (face*)malloc(sizeof(face));
        data->point_head->next = NULL;
        data->face_head->next = NULL;
}

void deallocate_filedata(filedata *data)
{
        deallocate_point(data->point_head);
        deallocate_face(data->face_head);
        data->point_head = NULL;
        data->face_head = NULL;
}

void deallocate_point(_point *node)
{
        if (node->next != NULL) {
                deallocate_point(node->next);
        } 
        free(node);
}

void deallocate_face(face *node)
{
        if (node->next != NULL) {
                deallocate_face(node->next);
        } 
        free(node);
}

void print_filedata(filedata *data) 
{
        _point *point_node = data->point_head->next;
        face *face_node = data->face_head->next;       
        while (point_node != NULL) {
                printf("V(%f, %f, %f)\n", point_node->x,
                                point_node->y, point_node->z);
                point_node = point_node->next;
        }
        while (face_node != NULL) {
                printf("F(%i, %i, %i)\n", face_node->v[0],
                                face_node->v[1], face_node->v[2]);
                face_node = face_node->next;
        }
        printf("%i vertices, %i faces\n", data->number_of_points,
                        data->number_of_triangles);
}
