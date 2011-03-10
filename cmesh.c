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
        sort_cyclic_order_at_vertices(m);
}

/* Free memory allocated in a previously initialized mesh. */
void deallocate_mesh(mesh *m) {
        free(m->vertices);
        free(m->edges);
        free(m->triangles);
        free(m->coordinates);
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

/* Sorts the incidence data at each vertex to respect the inherited
 * cyclic ordering. In other words, if two edges e_j and e_k 
 * (or vertices or triangles) are listed in position i and i+1, then  
 * e_k is the next edge encountered when rotating counter-clockwise from
 * edge e_j. For vertices on the bounday, the first and last edges listed
 * are boundary edges.
 */
void sort_cyclic_order_at_vertices(mesh *m)
{
        vertex *v;
        void *iv[MAX_DEGREE * 2 + 1];
        void *ie[MAX_DEGREE * 2 + 1];
        int head, tail, i;
        
        for (i=0; i < m->rank[0]; i++) {
                v = &(m->vertices[i]);
                head = MAX_DEGREE;
                tail = MAX_DEGREE;
                ie[head] = v->incident_edges[0];
                iv[head] = v->incident_vertices[0];
                while (tail - head < v->degree-1) {
                      head -= find_clockwise_edge(v, 
                                      &(ie[head]), &(iv[head]));
                      tail += find_counterclockwise_edge(v, 
                                      &(ie[tail]), &(iv[tail]));
                }
                sort_incident_edges_and_vertices(v, 
                                &(ie[head]), &(iv[head]));
                sort_incident_triangles(v);
        }
}

/* Finds the edge incident to v, immediately clockwise from ie[0],
 * and records it in ie[-1]. Does the same for the vertices.
 * If there is no edge clockwise from ie[0] (it is on the boundary)
 * then the function returns 0. Otherwise, it returns 1.
 */
int find_clockwise_edge(vertex *v, void *ie[], void *iv[])
{
        int i;
        edge *e = (edge *)ie[0];
        triangle *t;
        if (e->vertices[0] == v) {
                t = (triangle *)(e->incident_triangles[0]);
        } else {
                t = (triangle *)(e->incident_triangles[1]);
        }
        if (!t) {
                return 0;
        }
        ie--;
        iv--;
        for (i=0; i<3; i++) {
                if (t->vertices[i] == v) {
                        (*ie) = t->edges[(i+2) % 3];
                        (*iv) = t->vertices[(i+1) % 3];
                        return 1;
                }          
        }
        printf("Error in find_clockwise_edge\n");
        exit(1);
}

/* Finds the edge incident to v, immediately counterclockwise from ie[0],
 * and records it in ie[1]. Does the same for the vertices.
 * If there is no edge counterclockwise from ie[0] (it is on the boundary)
 * then the function returns 0. Otherwise, it returns 1.
 */
int find_counterclockwise_edge(vertex *v, void *ie[], void *iv[])
{
        int i;
        edge *e = (edge *)ie[0];
        triangle *t;
        if (e->vertices[0] == v) {
                t = (triangle *)(e->incident_triangles[1]);
        } else {
                t = (triangle *)(e->incident_triangles[0]);
        }
        if (!t) {
                return 0;
        }
        ie++;
        iv++;
        for (i=0; i<3; i++) {
                if (t->vertices[i] == v) {
                        (*ie) = t->edges[(i+1) % 3];
                        (*iv) = t->vertices[(i+2) % 3];
                        return 1;
                }          
        }
        printf("Error in find_counterclockwise_edge\n");
        exit(1);
}

/* Copies list of pointers ie to incident_edges and iv to 
 * incident_vertices. The assumption is that ie and iv are sorted
 * lists which just need to be recorded.
 */
void sort_incident_edges_and_vertices(vertex *v, void *ie[], void *iv[])
{
        int i;
        for (i=0; i<v->degree; i++) {
                v->incident_edges[i] = ie[i];
                v->incident_vertices[i] = iv[i];
        }
}

/* Assuming that the incident vertices and edges have been cyclically
 * ordered, this function does the same to the incident triangles.
 */
void sort_incident_triangles(vertex *v)
{
        int i;
        edge *e1, *e2;
        for (i=0; i < v->degree - v->boundary; i++) {
                e1 = (edge *)v->incident_edges[i];
                e2 = (edge *)v->incident_edges[(i+1) % v->degree];
                v->incident_triangles[i] = get_common_triangle(e1, e2);
        }
}

/* Returns a pointer to the triangle containing edges e1 and e2.
 * Exits (ungracefully) if the edges are not common to a triangle.
 */
void *get_common_triangle(edge *e1, edge *e2)
{
        int i, j;
        triangle *t;
        for (i=0; i<2; i++) {
                t = (triangle *)(e1->incident_triangles[i]);
                if (!t) {
                        for (j=0; j<3; j++) {
                                if (t->edges[j] == e2)
                                        return (void *)t;
                        }
                }
        } 
        printf("Edges do not share a common triangle.\n");
        exit(1);
}

/* returns the coordinate location of vertex v
 */
point* get_coordinate(mesh *m, vertex *v)
{
        int i = v->index;
        return &(m->coordinates[i]);
}

/* Calculates the Euclidean distance between two points in 3-space.
 */
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
 * to be incident. If so, a pointer to the edge joining them
 * is returned. If not, the NULL pointer is returned.
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

/* Adds index information to mesh data structs. 
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
