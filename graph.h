// graph.h

void calc_vertex_boundary_distances(mesh *m);
void create_active_list(mesh *m, int b, std::list<vertex*>& active);
void run_through_active_list(mesh *m, int b, std::list<vertex*>& active);
void calc_vertex_config(vertex* v, const geodesic &g, bool on_boundary=false);
void calc_next_vertex_geodesic(vertex *v, int i, geodesic &g, int b);
void add_geodesic_to_vertex(vertex *v, const geodesic &g);
