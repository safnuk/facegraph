// cmesh.h


#define max_degree  15

/* Struct which encodes incidence data at a vertex. The list of
 * incident vertices should be in alignment with the list of incident
 * edges. i.e. the first incident vertex should be connected to
 * the current vertex by the first incident edge, and so on.
 *
 * Also, all incidence data is cyclically ordered (counterclockwise).
 *
 * degree is the number of incident edges and vertices
 * degree - bounday is the number of incident triangles,
 * so boundary is 1 if vertex in on the boundary, 0 otherwise.
 *
 * dtheta_du records the derivatives of the inner angles of
 * all triangles incident to the vertex (s = e^u). They are pointers
 * to previously calculated Hessian data (stored in mesh struct).
 *
 * s = (e^r - 1) / (e^r + 1), where r is the radius from the circle
 * packing metric.
 */
typedef struct {
        int degree;    
        int boundary;
        void *incident_vertices[max_degree];
        void *incident_triangles[max_degree];
        void *incident_edges[max_degree];
        void *link_edges[max_degree];
        double *dtheta_du[max_degree][3];
        double s;
        double s0;
        int index;
} vertex;

/* Struct for encoding edge connectivity. Vertices are the indices to
 * the two endpoints of the edge. If vertices are listed in the order
 * (i, j), then triangle 1 should be on the left when traveling from 
 * left to right. i.e. it agrees with the ordering of triangle 1, 
 * and is opposite to the ordering of triangle 2 (if the edge is 
 * interior).
 *
 * cos_angle stores the cosine of the intersection of the two
 * circles (centered at the end points of the edge) used in
 * the circle packing metric.
 *
 * cosh_length stores the hyp cosine of the length of the edge,
 * computed using the circle packing metric.
 */
typedef struct {
        void *vertices[2];
        void *incident_triangles[2];
        double cos_angle;
        double cosh_length;
        double sinh_length;
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
 * means edges should be listed in the order (A, B, C).
 *
 * Inner angles are measured using hyperbolic triangles laws, 
 * with edge lengths coming from the circle packing metric.
 *
 * hessian is the matrix of derivatives [d theta_i / du_j].
 */
typedef struct {
        void *vertices[3];
        void *edges[3];
        double inner_angles[3];
        double cos_angles[3];
        double sin_angles[3];
        double hessian[3][3];
        int index;
} triangle;

/* Struct recording the (x, y, z) coordinates of a point.
 * This is mostly to avoid having to allocate 2D arrays.
 */
typedef struct {
        double x;
        double y;
        double z;
} point;

/* Struct encoding incidence data of a triangulated surface.
 * All cross references (eg incident_vertices in a vertex) are
 * indices to the arrays of objects.
 *
 * ranks is the numbers of vertices, edges and triangles.
 * coordinates is an array of (x,y,z) triples giving the
 * position of each vertex.
 *
 * The quantity f represents the integral of Sum(K_i du_i)
 * from the basepoint given by the initial circle packing
 * metric to the current state of the circle packing metric.
 * The gradient of f is the curvatures of the mesh, so 
 * minimizing f is equivalent to finding 0 curvature.
 *
 * TODO: Check triangulation data for
 *      - isolated vertices
 *      - bivalent vertices
 *      - improvable trivalent and quadvalent vertices
 */
typedef struct {
        double f;
        int ranks[3];  // ranks of i-th chain groups
        vertex *vertices;
        edge *edges;
        triangle *triangles;
        point *coordinates;
} mesh;

int initialize_mesh(mesh *m, filedata *fd);
void deallocate_mesh(mesh *m);
void add_indices(mesh *m);
void copy_points(_point *head, mesh *m);
void construct_triangles(face *head, mesh *m);
int construct_edges(mesh *m);
void *get_incident_edge(vertex *v1, vertex *v2);
void add_incident_vertices_and_edge(vertex *v1, vertex *v2, edge *e);
int valid_pointers(mesh *m);
void sort_cyclic_order_at_vertices(mesh *m);
int find_clockwise_edge(vertex *v, void *ie[], void *iv[]);
int find_counterclockwise_edge(vertex *v, void *ie[], void *iv[]);
void sort_incident_edges_and_vertices(vertex *v, void *ie[], void *iv[]);
void sort_incident_triangles(vertex *v);
void *get_common_triangle(edge *e1, edge *e2);
void add_link_edges(mesh *m);
void construct_vertex_hessian_pointers(mesh *m);
int get_vertex_position_in_triangle(vertex *v, triangle *t);
point* get_coordinate(mesh *m, vertex *v);

double calc_distance(point *p1, point *p2);
double update_f_and_s(mesh *m, double *u, double ds);
void calc_curvatures(mesh *m, double *K);
void calc_inner_angles(mesh *m);
double calc_curvature(vertex *v);
double curvature_integrand(double s, void *instance);
void calc_edge_length (edge *e);
double min(double *x, int n);
double max(double *x, int n);

void calc_hessian(mesh *m);
void calc_dtheta_dl(triangle *t, double A[3][3]);
void calc_dl_ds(triangle *t, double A[3][3]);
void calc_ds_du(triangle *t, double A[3][3]);
void calc_3_matrix_product(double A[3][3], double B[3][3],
    double C[3][3], double D[3][3]);

void print_mesh(mesh *m);
void print_coordinate(point *p);
void print_vertex(vertex *v);
void print_edge(edge *e);
void print_triangle(triangle *t);
