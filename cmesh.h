// cmesh.h


#define max_degree  15
#define max_boundaries 4

struct edge;
struct triangle;

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
struct vertex {
  public:
        vertex() : geodesics() {}
        int degree;    
        int boundary;
        vertex* incident_vertices[max_degree];
        triangle* incident_triangles[max_degree];
        edge* incident_edges[max_degree];
        edge* link_edges[max_degree];
        double* dtheta_du[max_degree][3];
        double* inner_angles[max_degree];
        double s;
        double s0;
        int index;
        std::list<geodesic> geodesics;
        geodesic shortest_paths[max_boundaries];
        geodesic shortest_path;
        vertex_config vc;
};


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
struct edge {
        vertex* vertices[2];
        triangle* incident_triangles[2];
        double cos_angle;
        double cosh_length;
        double sinh_length;
        int index;
};

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
struct triangle {
        vertex* vertices[3];
        edge* edges[3];
        double inner_angles[3];
        double cos_angles[3];
        double sin_angles[3];
        double hessian[3][3];
        int index;
};

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
        filedata* fd;
        double f;
        int ranks[3];  // ranks of i-th chain groups
        int boundary_count;
        vertex *vertices;
        edge *edges;
        triangle *triangles;
        point *coordinates;
        int* boundary_edges;
        std::list<edge*> boundary_cycles[max_boundaries];
        double boundary_lengths[max_boundaries];
        vertex_config vc;
        CompRow_Mat_double hessian;
} mesh;

int initialize_mesh(mesh *m, filedata *fd);
void construct_simplices(mesh *m, filedata* fd);
void fix_bivalent_vertices(mesh* m, filedata* fd);
void remove_isolated_vertices(mesh* m, filedata* fd);
void remap_triangle_vertices(filedata* fd, int* vertex_mapping);
triangle* get_other_incident_triangle(edge* e, triangle* t);
void list_triangle_vertices(triangle* t, int* vertices, int offset);
void construct_edge_bisectors(mesh* m, filedata* fd,
                              triangle* t, int* vertices, int offset);
face* add_new_triangles(face* face_node, int* t1, int* t2, int* t3);
void construct_new_triangle_list(filedata* fd, face* face_head,
                                 int* triangles_to_keep);
void double_mesh(mesh *m, mesh *m_double);
void split_doubled_mesh(mesh *m, mesh *m_double);
void deallocate_mesh(mesh *m);
void add_indices(mesh *m);
void copy_points(_point *head, mesh *m);
void construct_triangles(face *head, mesh *m);
int construct_edges(mesh *m);
void calc_boundaries(mesh *m);
int get_next_component_start(int* track_boundary_edges, int n);
edge* get_next_boundary_edge(edge* e);
void *get_incident_edge(vertex *v1, vertex *v2);
void add_incident_vertices_and_edge(vertex *v1, vertex *v2, edge *e);
int valid_pointers(mesh *m);
void double_vertices(mesh *m, mesh *m_double);
void double_edges(mesh *m, mesh *m_double);
void double_triangles(mesh *m, mesh *m_double);
void sort_cyclic_order_at_vertices(mesh *m);
int find_clockwise_edge(vertex* v, edge* ie[], vertex* iv[]);
int find_counterclockwise_edge(vertex* v, edge* ie[], vertex* iv[]);
void sort_incident_edges_and_vertices(vertex* v, edge* ie[], vertex* iv[]);
void sort_incident_triangles(vertex *v);
triangle* get_common_triangle(edge *e1, edge *e2);
void add_link_edges(mesh *m);
void construct_vertex_hessian_pointers(mesh *m);
void construct_vertex_inner_angle_pointers(mesh *m);
int get_vertex_position_in_triangle(vertex *v, triangle *t);
int get_edge_position_at_vertex(edge* e, vertex* v);
bool edges_share_triangle(edge* e1, edge* e2);
point* get_coordinate(mesh *m, vertex *v);

double calc_distance(point *p1, point *p2);
double update_f_and_s(mesh *m, double *s, int n);
void update_s_and_edge_lengths(mesh *m, MV_Vector_double &s);
void calc_curvatures(mesh *m, MV_Vector_double &K);
void calc_inner_angles(mesh *m);
double calc_curvature(vertex *v);
double curvature_integrand(double s, void *instance);
void calc_edge_lengths(mesh *m);
void calc_edge_length (edge *e);
double min(double *x, int n);
double max(double *x, int n);
void calc_boundary_lengths(mesh* m);

void clear_geodesic_lists(mesh* m);

void calc_hessian(mesh *m);
void calc_dtheta_dl(triangle *t, double A[3][3]);
void calc_dl_ds(triangle *t, double A[3][3]);
void calc_ds_du(triangle *t, double A[3][3]);
void calc_3_matrix_product(double A[3][3], double B[3][3],
    double C[3][3], double D[3][3]);
int calc_number_of_nonzero_hessian_entries(mesh *m);

void print_mesh(mesh *m);
void print_coordinate(point *p);
void print_vertex(vertex *v);
void print_edge(edge *e);
void print_triangle(triangle *t);

bool is_number(double x);
