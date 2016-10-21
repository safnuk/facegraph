// mesh_improver.h
#ifndef MeshImprover
#define MeshImprover

int find_problem_vertices(mesh* m, filedata* fd);
void list_incident_vertices(vertex* v, int* incident_points, 
                            int* vertex_hit_list);
void remove_triangles_from_keep_list(vertex* v, int* triangles_to_keep);
void construct_incident_edge_bisectors(mesh* m, filedata* fd, vertex* v, 
                                       int* new_middle_points);
face* add_new_triangles_for_trivalent(face* face_node, int p, int* v, int* m,
                                      int degree, int boundary);
#endif
