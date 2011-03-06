// cmesh.h


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

int initialize_mesh(mesh *m, filedata *fd);
void deallocate_mesh(mesh *m);
void add_indices(mesh *m);
void print_mesh(mesh *m);
void print_vertex(vertex *v);
void print_edge(edge *e);
void print_triangle(triangle *t);
void copy_points(_point *head, mesh *m);
void construct_triangles(face *head, mesh *m);
int construct_edges(mesh *m);
void *get_incident_edge(vertex *v1, vertex *v2);
void add_incident_vertices_and_edge(vertex *v1, vertex *v2, edge *e);
int valid_pointers(mesh *m);
