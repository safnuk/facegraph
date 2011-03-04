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
typedef struct _vertex {
        int degree;    
        int boundary;
        struct _vertex incident_vertices[max_degree];
        int incident_triangles[max_degree];
        int incident_edges[max_degree];
} vertex;

/* Struct for encoding edge connectivity. Vertices are the indices to
 * the two endpoints of the edge. If vertices are listed in the order
 * (i, j), then triangle 1 should be on the left when traveling from 
 * left to right. i.e. it agrees with the ordering of triangle 1, 
 * and is opposite to the ordering of triangle 2 (if the edge is 
 * interior).
 */
typedef struct {
        int vertices[2];
        int incident_triangles[2];
} edge;

int get_other_vertex(edge *e, int v);

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
        int vertices[3];
        int edges[3];
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
void copy_points(float points[][3], mesh *m);
void construct_triangles(int faces[][3], mesh *m);
int construct_edges(mesh *m);
void get_vertex_pair(mesh *m, triangle *t, int j, vertex *v1, vertex *v2,
                int &vi1, int &vi2);
int valid_pointers(mesh *m);

int main()
{
        int i;
        vertex v;
        mesh m;
        float points[10][3];
        int faces[15][3];
        v.degree = 2;
        v.boundary = 0;
        v.incident_vertices[0] = 10;
        v.incident_vertices[1] = 15;
        for(i=0; i<v.degree; i++) {
                printf("%i\n", v.incident_vertices[i]);
        }
        initialize_mesh(&m, points, faces, 10, 15);

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
                        t->vertices[j] = faces[i][j];
                        // Add incidence data to vertex
                        v = m->vertices[faces[i][j]];
                        ti = v->degree - v->boundary;
                        v->incident_triangles[ti] = i;
                        (v->boundary)--;
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
        int ti, vi1, vi2, j, k;
        int duplicate_edge;
        triangle *t;
        edge *e;
        vertex *v1, *v2;

        for (ti=0; ti< m->ranks[2]; ti++) {
                t = &(m->triangles[ti]);
                for (j=0; j < 3; j++) {
                        get_vertex_pair(m, t, j, v1, v2, vi1, vi2);
                        duplicate_edge = 0;
                        for (k=0; k < v2->degree; k++) {

                        }
                }
        }
        return edge_count;
}

void get_vertex_pair(mesh *m, triangle *t, int j, vertex *v1, vertex *v2,
                int &vi1, int &vi2)
{
        vi1 = t->vertices[j];
        vi2 = t->vertices[(j+1) % 3];
        v1 = &(m->vertices[vi1]);
        v2 = &(m->vertices[vi2])
}

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
int get_other_vertex(edge *e, int v)
{
        if (v == e->vertices[0])
                return e->vertices[1];
        if (v == e->vertices[1])
                return e->vertices[0];
        else 
                return -1;
}
