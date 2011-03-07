// cmesh.c

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cfile_io.h"
#include "cmesh.h"



/* Given a filedata structure (read from a .mesh file), 
 * the function calculates all incidence 
 * relations for the mesh data struct.
 *
 * Point triples in a face should be listed in counter-clockwise
 * order, and are indices to the inherited order from the points
 * array.
 */
int initialize_mesh(mesh *m, filedata *fd)
{
        int max_number_of_edges = fd->number_of_points +
                fd->number_of_triangles + 10; 
        int edge_count = 0;
        m->ranks[0] = fd->number_of_points;
        m->ranks[2] = fd->number_of_triangles;
        m->vertices = (vertex*) malloc(m->ranks[0] * sizeof(vertex));
        m->edges = (edge*) malloc(max_number_of_edges * sizeof(edge));
        m->triangles = (triangle*) malloc(m->ranks[2] * sizeof(triangle));
        m->coordinates = (point*)malloc(m->ranks[0] * sizeof(point));
        m->edge_lengths = (double*)malloc(max_number_of_edges *
                        sizeof(double)); 
        m->circle_angles = (double*)malloc(max_number_of_edges *
                        sizeof(double)); 
        m->radii = (double*)malloc(m->ranks[0] * sizeof(double));
        if (!valid_pointers(m)) {
                printf("Memory allocation failure!\n Goodbye.\n");
                deallocate_mesh(m);
                exit;
        }
        copy_points(fd->point_head, m);
        construct_triangles(fd->face_head, m);
        edge_count = construct_edges(m);
        m->ranks[1] = edge_count;
        add_indices(m);
}

/* Free memory allocated in a previously initialized mesh. */
void deallocate_mesh(mesh *m) {
        free(m->vertices);
        free(m->edges);
        free(m->triangles);
        free(m->coordinates);
        free(m->edge_lengths);
        free(m->circle_angles);
        free(m->radii);
}
/* Copy array of (x,y,z) coordinate data into mesh's
 * coordinates array.
 */
void copy_points(_point *head, mesh *m)
{
        _point *node = head;
        int i = 0;
        while (node->next != NULL) {
                node = node->next;
                m->coordinates[i].x = node->x;
                m->coordinates[i].y = node->y;
                m->coordinates[i].z = node->z;
                i++;
        }
}


/* Initializes triangle data for mesh, and adds incidence data 
 * to vertices.
 */
void construct_triangles(face *head, mesh *m)
{
        int i=0; 
        int j, ti;
        triangle *t;
        vertex *v;
        face *node = head;

        while (node->next != NULL) {
                node = node->next;
                t = &(m->triangles[i]);
                for (j=0; j<3; j++) {
                        // Add incidence data to vertex and triangle
                        v = &(m->vertices[node->v[j]]);
                        ti = v->degree - v->boundary;
                        v->incident_triangles[ti] = (void *)t;
                        (v->boundary)--;
                        t->vertices[j] = (void *)v;
                }
                i++;
        }
}

/* Initializes edge data for mesh, and adds incidence data
 * to vertices and triangles.
 * Returns the true count of number of edges.
 */
int construct_edges(mesh *m)
{
        int edge_count = 0;
        int i, j, opposite_edge, k;
        triangle *t;
        edge *e;
        vertex *v1, *v2;

        for (i=0; i< m->ranks[2]; i++) {
                t = &(m->triangles[i]);
                for (j=0; j < 3; j++) {
                        v1 = (vertex *)(t->vertices[j]);
                        v2 = (vertex *)(t->vertices[(j+1) % 3]);
                        e = (edge *)get_incident_edge(v1, v2);
                        if (e != NULL) {
                                e->incident_triangles[1] = (void *)t;
                        }
                        else {
                                e = &(m->edges[edge_count]);
                                edge_count++;
                                e->incident_triangles[0] = (void *)t;
                                e->incident_triangles[1] = NULL;
                                add_incident_vertices_and_edge(v1, v2, e);
                        }
                        opposite_edge = (j+2) % 3;
                        t->edges[opposite_edge] = (void *)e;
                }
        }
        return edge_count;
}

/* returns the coordinate location of vertex v
 */
point* get_coordinate(mesh *m, vertex *v)
{
        int i = v->index;
        return &(m->coordinates[i]);
}

double calc_distance(point *p1, point *p2)
{
        double diff;
        double sum;
        diff = p1->x - p2->x;
        sum = diff * diff;
        diff = p1->y - p2->y;
        sum += diff * diff;
        diff = p1->z - p2->z;
        sum += diff * diff;
        return sqrt(sum);
}


// TODO: Update comment.
/* Checks to see if vertices v1 and v2 have been previously found
 * to be incident. If so, edge e is set to the edge joining them.
 */
void * get_incident_edge(vertex *v1, vertex *v2) 
{
        int i;
        vertex *v;
        for (i=0; i < v1->degree; i++) {
                v = (vertex *)(v1->incident_vertices[i]);
                if (v == v2) {
                        return v1->incident_edges[i];
                } 
        }
        return NULL;
}

/* Adds incidence relation of v1 to v2 and vice versa. 
 * Assumes that edge e joins v1 to v2, so adds its incidence
 * data to v1 and v2 as well.
 */
void add_incident_vertices_and_edge(vertex *v1, vertex *v2, edge *e)
{
        v1->incident_vertices[v1->degree] = (void *)v2;
        v1->incident_edges[v1->degree] = (void *)e;
        (v1->degree)++;;
        (v1->boundary)++;
        v2->incident_vertices[v2->degree] = (void *)v1;
        v2->incident_edges[v2->degree] = (void *)e;
        (v2->degree)++;;
        (v2->boundary)++;
        e->vertices[0] = (void *)v1;
        e->vertices[1] = (void *)v2;
}

/* Checks to see if all arrays allocated in initialize_mesh are valid pointers. */
int valid_pointers(mesh *m)
{
        if (m->vertices == NULL)
                return 0;
        if (m->edges == NULL)
                return 0;
        if (m->triangles == NULL)
                return 0;
        if (m->coordinates == NULL)
                return 0;
        if (m->edge_lengths == NULL)
                return 0;
        if (m->radii == NULL)
                return 0;
        if (m->circle_angles == NULL)
                return 0;
        return 1;
}

/* Finds the vertex different from v, incident to edge e.
 * If v is not incident to e, then -1 is returned.
 */
/*
int get_other_vertex(edge *e, int v)
{
        if (v == e->vertices[0])
                return e->vertices[1];
        if (v == e->vertices[1])
                return e->vertices[0];
        else 
                return -1;
}
*/

/* MAY NOT BE NEEDED.
 * Adds index information to mesh data structs. 
 * i.e. records position of each object within its
 * respective array.
 */
void add_indices(mesh *m) 
{
        int i;
        vertex *v;
        edge *e;
        triangle *t;
        for (i=0; i < m->ranks[0]; i++) {
                v = &(m->vertices[i]);
                v->index = i;
        }
        for (i=0; i < m->ranks[1]; i++) {
                e = &(m->edges[i]);
                e->index = i;
        }
        for (i=0; i < m->ranks[2]; i++) {
                t = &(m->triangles[i]);
                t->index = i;
        }
}

void print_mesh(mesh *m) 
{
        int i;

        printf("==== Coordinates ====\n");
        for (i=0; i < m->ranks[0]; i++) {
                printf("[%i]", i);
                print_coordinate(&(m->coordinates[i]));
                printf("\n");
        }
        printf("==== Vertices ====\n");
        for (i=0; i < m->ranks[0]; i++) {
                 print_vertex(&(m->vertices[i]));
        }
        printf("==== Edges ====\n");
        for (i=0; i < m->ranks[1]; i++) {
                print_edge(&(m->edges[i]));
        }
        printf("==== Triangles ====\n");
        for (i=0; i < m->ranks[2]; i++) {
                print_triangle(&(m->triangles[i]));
        }
        printf("Ranks = (%i, %i, %i)", m->ranks[0],
                        m->ranks[1], m->ranks[2]);
        printf("  Euler = %i\n", m->ranks[0] - 
                        m->ranks[1] + m->ranks[2]);
}

void print_coordinate(point *p)
{
        printf("(%f, %f, %f)  ", p->x, p->y, p->z);
}

void print_vertex(vertex *v) 
{
        int j, k;
        vertex *v1;
        edge *e;
        triangle *t;

        printf("(V%i)  V:[", v->index);
        for (j=0; j < v->degree; j++) {
                v1 = (vertex *)(v->incident_vertices[j]);
                k = v1->index;
                printf(" %i ", k);
        }
        printf("]  E:[");
        for (j=0; j < v->degree; j++) {
                e = (edge *)(v->incident_edges[j]);
                k = e->index;
                printf(" %i ", k);
        }
        printf("]  T:[");
        for (j=0; j < v->degree - v->boundary; j++) {
                t = (triangle *)(v->incident_triangles[j]);
                k= t->index;
                printf(" %i ", k);
        }
        printf("] Deg:%i, Boundary:%i\n", v->degree, v->boundary);
}

void print_edge(edge *e) 
{
        int v1, v2;
        int t1, t2;
        void *t;
        v1 = ((vertex *)(e->vertices[0]))->index;
        v2 = ((vertex *)(e->vertices[1]))->index;
        t1 = ((triangle *)(e->incident_triangles[0]))->index;
        t = e->incident_triangles[1];
        printf("(E%i) ", e->index);
        if (t == NULL) {
                printf("V:[%i, %i]  T:[%i]\n", v1, v2, t1);
        }
        else {
                t2 = ((triangle *)t)->index;
                printf("V:[%i, %i]  T:[%i, %i]\n", v1, v2, t1, t2);
        }
}

void print_triangle(triangle *t)
{
        int v[3];
        int e[3];
        int i;

        for(i=0; i < 3; i++) {
                v[i] = ((vertex *)(t->vertices[i]))->index;
                e[i] = ((edge *)(t->edges[i]))->index;
        }
        printf("(T%i)  V:[%i, %i, %i]", t->index, v[0], v[1], v[2]);
        printf(" E:[%i, %i, %i]\n", e[0], e[1], e[2]);
}
