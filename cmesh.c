// cmesh.c

#include <stdlib.h>
#include <stdio.h>

#define max_degree  15

/* Struct which encodes incidence data at a vertex. The list of
 * incident vertices should be in alignment with the list of incident
 * edges. i.e. the first incident vertex should be connected to
 * the current vertex by the first incident edge, and so on.
 *
 * degree is the number of incident edges and vertices
 * degree - bounday is the number of incident triangles,
 * so boundary is 1 if vertex in on the boundary, 0 otherwise.
 *
 * TODO: Make lists respect inherent cyclic ordering.
 */
typedef struct {
        int degree;    
        int boundary;
        void *incident_vertices[max_degree];
        void *incident_triangles[max_degree];
        void *incident_edges[max_degree];
        int index;
} vertex;

/* Struct for encoding edge connectivity. Vertices are the indices to
 * the two endpoints of the edge. If vertices are listed in the order
 * (i, j), then triangle 1 should be on the left when traveling from 
 * left to right. i.e. it agrees with the ordering of triangle 1, 
 * and is opposite to the ordering of triangle 2 (if the edge is 
 * interior).
 */
typedef struct {
        void *vertices[2];
        void *incident_triangles[2];
        int index;
} edge;

// int get_other_vertex(edge *e, vertex *v);

/* Triangle struct used to encode connectivity structure.
 * Vertices should be a counter-clockwise list of indices.
 * Edge indices are also listed counter-clockwise, with the
 * edge being opposite to the vertex.
 *
 * i.e. if a triangle has vertices a, b, c with opposite 
 * edges A, B, C, then vertices listed in the order (a, b, c) 
 * means edges should be listed in the order (A, B, C)
 */
typedef struct {
        void *vertices[3];
        void *edges[3];
        int index;
} triangle;

/* Struct recording the (x, y, z) coordinates of a point.
 * This is mostly to avoid having to allocate 2D arrays.
 */
typedef struct {
        float x;
        float y;
        float z;
} point;

/* Struct encoding incidence data of a triangulated surface.
 * All cross references (eg incident_vertices in a vertex) are
 * indices to the arrays of objects.
 *
 * ranks is the numbers of vertices, edges and triangles.
 * coordinates is an array of (x,y,z) triples giving the
 * position of each vertex.
 *
 * TODO: Check triangulation data for
 *      - isolated vertices
 *      - bivalent vertices
 *      - improvable trivalent and quadvalent vertices
 */
typedef struct {
        int ranks[3];  // ranks of i-th chain groups
        vertex *vertices;
        edge *edges;
        triangle *triangles;
        point *coordinates;
} mesh;

int initialize_mesh(mesh *m, float points[][3], int faces[][3],
                int number_of_points, int number_of_triangles);
void deallocate_mesh(mesh *m);
void add_indices(mesh *m);
void print_mesh(mesh *m);
void print_vertex(vertex *v);
void print_edge(edge *e);
void print_triangle(triangle *t);
void copy_points(float points[][3], mesh *m);
void construct_triangles(int faces[][3], mesh *m);
int construct_edges(mesh *m);
void * get_incident_edge(vertex *v1, vertex *v2);
void add_incident_vertices_and_edge(vertex *v1, vertex *v2, edge *e);
int valid_pointers(mesh *m);

int main()
{
        int i;
        vertex v;
        mesh m;
        float points[7][3] ={ 1,1,1, 1,1,1, 1,1,1,
              1,1,1, 1,1,1, 1,1,1, 1,1,1};
        int faces[6][3] = { 0,3,1, 1,3,2, 2,3,4, 2,4,5,
            0,6,3, 3,6,4};
        initialize_mesh(&m, points, faces, 7,  6);
        print_mesh(&m);
        deallocate_mesh(&m);
}

/* Given an array of point coordinates, and triples of points
 * which form triangles, the function calculates all incidence 
 * relations for the mesh data struct.
 *
 * Point triples in a face should be listed in counter-clockwise
 * order, and are indices to the inherited order from the points
 * array.
 */
int initialize_mesh(mesh *m, float points[][3], int faces[][3],
                int number_of_points, int number_of_triangles)
{
        int max_number_of_edges = number_of_points + number_of_triangles; 
        int edge_count = 0;
        m->ranks[0] = number_of_points;
        m->ranks[2] = number_of_triangles;
        m->vertices = (vertex*) malloc(number_of_points * sizeof(vertex));
        m->edges = (edge*) malloc(max_number_of_edges * sizeof(edge));
        m->triangles = (triangle*) malloc(number_of_triangles * 
                        sizeof(triangle));
        m->coordinates = malloc(number_of_points * sizeof(point));
        if (!valid_pointers(m)) {
                printf("Memory allocation failure!\n Goodbye.\n");
                deallocate_mesh(m);
                exit;
        }
        copy_points(points, m);
        construct_triangles(faces, m);
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
}
/* Copy array of (x,y,z) coordinate data into mesh's
 * coordinates array.
 */
void copy_points(float points[][3], mesh *m)
{
        int i;
        for (i=0; i < m->ranks[0]; i++) {
                m->coordinates[i].x = points[i][0];
                m->coordinates[i].y = points[i][1];
                m->coordinates[i].z = points[i][2];
        }
}


/* Initializes triangle data for mesh, and adds incidence data 
 * to vertices.
 */
void construct_triangles(int faces[][3], mesh *m)
{
        int i, j, ti;
        triangle *t;
        vertex *v;
        for (i=0; i < m->ranks[2]; i++) {
                t = &(m->triangles[i]);
                for (j=0; j<3; j++) {
                        // Add incidence data to vertex and triangle
                        v = &(m->vertices[faces[i][j]]);
                        ti = v->degree - v->boundary;
                        v->incident_triangles[ti] = (void *)t;
                        (v->boundary)--;
                        t->vertices[j] = (void *)v;
                }
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
