// graph.h
const int kVertexActive = 1;
const int kVertexNotActive = 0;
const double kErrorThreshold = 1e-6;

/* Struct used to store information relevant
 * for vertices of the mesh which are adjacent to
 * the cut-locus graph (graph of points equidistant
 * to two or more boundaries). More precisely, a
 * graph_vertex is an endpoint of a "transition edge",
 * i.e. an edge whose two vertices are closest to different
 * boundaries.
 *      v points to the actual vertex
 *      e points to the transition edge which gives rise to v
 *      edge_index is the index of e in v->incident_edges
 *      v_opposiste points to the other vertex in e
 */
struct graph_vertex {
  graph_vertex(vertex* v, vertex* v_opposite, edge* e, int edge_index) :
    v(v), v_opposite(v_opposite), e(e), edge_index(edge_index) {}
  vertex* v;
  vertex* v_opposite;
  edge* e;
  int edge_index;
};

struct ribbon_graph {
  int half_edge_count;
  std::vector<int> edge_perm;
  std::vector<int> face_perm;
  std::vector<double> metric;
};

void calc_cutlocus_graph(mesh* m, ribbon_graph* gamma);
void find_transition_edges(mesh* m, std::list<edge*>& transition_edges);
void find_and_sort_transition_vertices(mesh* m, std::list<edge*> const& transition_edges,
                std::list<graph_vertex>* transition_vertices);
void sort_transition_vertices(std::list<graph_vertex>& vertices);
void calc_graph_cycles(mesh* m, std::list<graph_vertex> const* transition_vertices, 
                ribbon_graph* gamma);
void calc_vertex_boundary_distances(mesh *m);
void create_active_list(mesh* m, int b, std::list<vertex*>& active);
void run_through_active_list(mesh* m, int b, std::list<vertex*>& active);
void calc_closest_boundaries(mesh* m);
int recalc_vertex_geodesic(vertex* v);
void calc_vertex_config(vertex* v, geodesic const& g, bool on_boundary);
void calc_next_vertex_geodesic(mesh* m, vertex* v, int k, geodesic& g, int b,
    geodesic const* g_source);
double calc_next_geodesic_edge_angle(vertex *v, vertex* v_next, double beta);
double calc_geodesic_edge_angle(vertex *v, const geodesic& path, int k);
void add_geodesic_to_vertex(vertex* v, vertex* orig_v, const geodesic& g);
geodesic calc_average_geodesic(vertex* v);
double normalize_angle(double angle);
double normalize_position(mesh* m, double position, int boundary);
bool compare_graph_vertices(graph_vertex const& gv1, graph_vertex const& gv2);

void partition_boundaries(int boundary_count, 
                          std::list<graph_vertex> const* transition_vertices,
                          std::list<double>* boundary_partitions);
void calc_half_edges_and_metric(mesh* m,
                                std::list<double> const* boundary_partitions,
                                ribbon_graph* gamma);
void calc_boundary_permutation(mesh* m,
                               std::list<double> const* boundary_partitions,
                               ribbon_graph* gamma);
void calc_edge_permutation(mesh const* m,
                           std::list<graph_vertex> const* transition_vertices,
                           std::list<double> const* boundary_partitions,
                           ribbon_graph* gamma);
int calc_half_edge_index(int boundary, double position, 
    std::list<double> const* boundary_partitions);

void print_ribbon_graph(ribbon_graph const* gamma);
void print_geodesic_list(vertex const* v);
