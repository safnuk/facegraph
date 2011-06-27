// graph.h

void calc_vertex_boundary_distances(mesh *m);
void create_active_list(mesh* m, int b, std::list<vertex*>& active);
void run_through_active_list(mesh* m, int b, std::list<vertex*>& active);
void calc_vertex_config(vertex* v, geodesic const& g, bool on_boundary=false);
void calc_next_vertex_geodesic(mesh* m, vertex* v, int k, geodesic& g, int b,
    geodesic const* g_source=NULL);
double calc_next_ geodesic_edge_angle(vertex *v, vertex* v_next, double beta);
double calc_geodesic_edge_angle(vertex *v, const geodesic& path, int k);
void add_geodesic_to_vertex(vertex* v, vertex* orig_v, const geodesic& g);
geodesic calc_average_geodesic(vertex* v);
double normalize_angle(double angle);
double normalize_position(mesh* m, double position, int boundary);
