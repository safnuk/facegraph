// cmesh.c

#define DEBUG_MODE 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>
#include <iostream>
#include <lbfgs.h>

#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"

#include "cfile_io.h"
#include "geodesic.h"
#include "csimpson.h"
#include "cmesh.h"
#include "mesh_improver.h"
#include "ccirclepack.h"

using namespace std;

const double edge_length_threshold = 1e-15;
const double radii_threshold = 1e-12;

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
  m->fd = fd;
  int max_number_of_edges = fd->number_of_points +
    fd->number_of_triangles + 10; 
  int edge_count = 0;
  m->ranks[0] = fd->number_of_points;
  m->ranks[2] = fd->number_of_triangles;
  m->vertices = new vertex[2* m->ranks[0]];
  m->edges = (edge*) malloc(2*max_number_of_edges * sizeof(edge));
  m->triangles = (triangle*) malloc(2*m->ranks[2] * sizeof(triangle));
  m->coordinates = (point*)malloc(2*m->ranks[0] * sizeof(point));
  if (!valid_pointers(m)) {
    printf("Memory allocation failure!\n Goodbye.\n");
    deallocate_mesh(m);
    exit;
  }
  construct_simplices(m, fd);
  fix_bivalent_vertices(m, fd);
  remove_isolated_vertices(m, fd);
  sort_cyclic_order_at_vertices(m);
  add_link_edges(m);
  m->f = 0;
  calc_boundaries(m);
  construct_vertex_hessian_pointers(m);
  construct_vertex_inner_angle_pointers(m);
}

void construct_simplices(mesh *m, filedata* fd)
{
  copy_points(fd->point_head, m);
  construct_triangles(fd->face_head, m);
  m->ranks[0] = fd->number_of_points;
  m->ranks[2] = fd->number_of_triangles;
  m->ranks[1] =  construct_edges(m);
  add_indices(m);
}

/* Looks for all triangles with bivalent vertices.
 * These triangles have their edges bisected and replaced
 * with three smaller triangles (the tip with the bivalent
 * vertex is removed). The adjacent triangle is cut 
 * into two smaller triangles.
 */
void fix_bivalent_vertices(mesh* m, filedata* fd)
{
  int old_boundary_triangle[3];
  int old_interior_triangle[3];
  int new_middle_points[3];
  int* triangles_to_keep = new int[m->ranks[2]];
  for (int i=0; i<m->ranks[2]; ++i) {
    triangles_to_keep[i] = 1;
  }
  face* face_head = (face*)malloc(sizeof(face));
  face_head->next = NULL;
  face* face_node = face_head;
  for (int i=0; i<m->ranks[0]; ++i) {
    if (m->vertices[i].degree == 2) {
      vertex* v = &(m->vertices[i]);
      triangle* t = v->incident_triangles[0];
      int k = get_vertex_position_in_triangle(v, t);
      triangle* boundary_triangle = t;
      triangle* interior_triangle = get_other_incident_triangle(t->edges[k], t);
      triangles_to_keep[boundary_triangle->index] = 0;
      triangles_to_keep[interior_triangle->index] = 0;
      list_triangle_vertices(boundary_triangle, old_boundary_triangle, ((k+1)%3));
      v = t->vertices[(k+1)%3];
      int k2 = get_vertex_position_in_triangle(v, interior_triangle);
      list_triangle_vertices(interior_triangle, old_interior_triangle, k2);
      construct_edge_bisectors(m, fd, t, new_middle_points, k);
      face_node = add_new_triangles(face_node, old_boundary_triangle, 
          old_interior_triangle, new_middle_points);
    }
  }
  construct_new_triangle_list(fd, face_head, triangles_to_keep);
  construct_simplices(m, fd);
  delete [] triangles_to_keep;

}

/* Assuming that triangle t is adjacent to edge e, function
 * returns the other triangle adjacent to e. If triangle t
 * is not adjacent to edge e, function fails badly.
 */
triangle* get_other_incident_triangle(edge* e, triangle* t)
{
  if (e->incident_triangles[0] == t) {
    return e->incident_triangles[1];
  }
  else if (e->incident_triangles[1] == t){
    return e->incident_triangles[0];
  } else {
    printf("Triangle not adjacent to edge in get_other_incident_triangle.\n");
    exit(1);
  }
}

/* Assuming that triangle t is comprised of vertices v0, v1, v2
 * (with indices i0, i1, i2), function stores the indices in array
 * vertices. The first vertex stored is specified by variable offset.
 */
void list_triangle_vertices(triangle* t, int* vertices, int offset)
{
  for (int i=0; i<3; ++i) {
    vertices[i] = t->vertices[(i+offset)%3]->index;
  }
}

/* The three edges of triangle t are bisected. Resulting midpoints are added
 * to the end of the list of vertices in fd, with the resulting indices stored
 * in array vertices. Variable offset is used to specify the position of the
 * bivalent vertex in triangle t.
 */
void construct_edge_bisectors(mesh* m, filedata* fd, triangle* t, int* vertices, int offset)
{
  _point* point_node = fd->point_head;
  while (point_node->next != NULL) {
    point_node = point_node->next;
  }
  for (int i=0; i<3; ++i) {
    edge* e = t->edges[(i+offset)%3];
    int j0 = e->vertices[0]->index;
    int j1 = e->vertices[1]->index;
    double x = (m->coordinates[j0].x + m->coordinates[j1].x) / 2;
    double y = (m->coordinates[j0].y + m->coordinates[j1].y) / 2;
    double z = (m->coordinates[j0].z + m->coordinates[j1].z) / 2;
    vertices[i] = fd->number_of_points;
    point_node = add_point_node(point_node, x, y, z);
    ++(fd->number_of_points);
  }
}

/* Assuming that a triangle with bivalent vertex is properly bisected,
 * functions constructs 5 new triangles that subdivide the bivalent triangle
 * and its adjacent triangle.
 */
face* add_new_triangles(face* face_node, int* t1, int* t2, int* t3)
{
  face_node = add_face_node(face_node, t1[0], t3[0], t3[2]);
  face_node = add_face_node(face_node, t3[0], t3[1], t3[2]);
  face_node = add_face_node(face_node, t3[0], t1[1], t3[1]);
  face_node = add_face_node(face_node, t2[0], t2[1], t3[0]);
  face_node = add_face_node(face_node, t3[0], t2[1], t2[2]);
  return face_node;
}

/* Function which constructs a new triangle list for fd by taking
 * all newly created triangles from the subdivision process and adding on all
 * the preexisitng triangles that were not subdivided.
 */
void construct_new_triangle_list(filedata* fd, face* face_head, int* triangles_to_keep)
{
  int count = 0;
  face* face_node = face_head;
  while(face_node->next != NULL) {
    ++count;
    face_node = face_node->next;
  }
  face* old_face_node = fd->face_head;
  bool delete_prior_node = true;
  for (int i=0; i<fd->number_of_triangles; ++i) {
    face* prior_face_node = old_face_node;
    old_face_node = old_face_node->next;
    if (delete_prior_node) {
      free(prior_face_node);
    }
    if (triangles_to_keep[i]) {
      delete_prior_node = false;
      face_node->next = old_face_node;
      face_node = face_node->next;
      ++count;
    } else {
      delete_prior_node = true;
    }
  }
  if (delete_prior_node) {
    free(old_face_node);
  }
  face_node->next = NULL;
  fd->number_of_triangles = count;
  fd->face_head = face_head;
}

/* Function which removes all the isolated vertices from
 * fd, and then reconstructs the mesh.
 */
void remove_isolated_vertices(mesh* m, filedata* fd)
{
  _point* point_node = fd->point_head;
  _point* prior_point_node = point_node;
  int* vertex_mapping = new int[m->ranks[0]];
  int keep_count = 0; 
  for (int i=0; i<m->ranks[0]; ++i) {
    point_node = prior_point_node->next;
    vertex* v = &(m->vertices[i]);
    if (v->degree == 0) {
      vertex_mapping[i] = -1;
      prior_point_node->next = point_node->next;
      free(point_node);
      point_node = prior_point_node;
    } else {
      vertex_mapping[i] = keep_count;
      ++keep_count;
      prior_point_node = point_node;
    }
  }
  fd->number_of_points = keep_count;
  remap_triangle_vertices(fd, vertex_mapping);
  construct_simplices(m, fd);
  delete [] vertex_mapping;
}

void remap_triangle_vertices(filedata* fd, int* vertex_mapping)
{
  face* face_node = fd->face_head->next;
  while (face_node != NULL) {
    for (int i=0; i<3; ++i) {
      face_node->v[i] = vertex_mapping[face_node->v[i]];
      if (face_node->v[i] == -1) {
        printf("Reference to isolated vertex in remap_triangle_vertices.\n");
        exit(1);
      }
    }
    face_node = face_node->next;
  }
}

void double_mesh(mesh *m, mesh *m_double)
{
  int i, j;
  int offset = 0;
  int* point_map = (int*) malloc(2 * m->ranks[0] * sizeof(int));
  filedata data;
  _point* point_node;
  face* face_node;
  point* p;
  int t[3];

  initialize_filedata(&data);
  point_node = data.point_head;
  face_node = data.face_head;
  if (!point_map) {
    printf("Memory allocation failure!\n Goodbye.\n");
    exit;
  }
  for (i=0; i < m->ranks[0]; i++) {
    point_node = add_point_node(point_node, m->coordinates[i].x,
        m->coordinates[i].y, m->coordinates[i].z);
    if (m->vertices[i].boundary) {
      point_map[i] = i;
    } else {
      point_map[i] = m->ranks[0] + offset;
      point_map[m->ranks[0] + offset] = i;
      offset++;
    }
  }
  for (i=m->ranks[0]; i < m->ranks[0] + offset; i++) {
    p = &(m->coordinates[point_map[i]]);
    point_node = add_point_node(point_node, p->x, p->y, p->z);
  }
  for (i=0; i < m->ranks[2]; i++) {
    for (j=0; j < 3; j++) {
      t[j] = ((vertex*)m->triangles[i].vertices[j])->index;
    }
    face_node = add_face_node(face_node, t[0], t[1], t[2]);
  }
  for (i=0; i < m->ranks[2]; i++) {
    for (j=0; j < 3; j++) {
      t[j] = point_map[((vertex*)m->triangles[i].vertices[j])->index];
    }
    face_node = add_face_node(face_node, t[0], t[2], t[1]);
  }
  data.number_of_points = m->ranks[0] + offset;
  data.number_of_triangles = 2 * m->ranks[2];

  initialize_mesh(m_double, &data);

  deallocate_filedata(&data);
  free(point_map);
}

void split_doubled_mesh(mesh *m, mesh *m_double)
{
  int i;

  for (i=0; i < m->ranks[0]; i++) {
    m->vertices[i].s = m_double->vertices[i].s;
  }
  for (i=0; i < m->ranks[1]; i++) {
    m->edges[i].cos_angle = m_double->edges[i].cos_angle;
    m->edges[i].cosh_length_minus1 = m_double->edges[i].cosh_length_minus1;
    m->edges[i].sinh_length = m_double->edges[i].sinh_length;
  }
}
/* Free memory allocated in a previously initialized mesh. */
void deallocate_mesh(mesh *m) {
  for (int i=0; i<m->boundary_count; ++i) {
    m->boundary_cycles[i].clear();
  }
  m->boundary_count=0;
  delete [] m->vertices;
  free(m->edges);
  free(m->triangles);
  free(m->coordinates);
  free(m->boundary_edges);
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
    m->vertices[i].degree = 0;
    m->vertices[i].boundary = 0;
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
      v->incident_triangles[ti] = t;
      (v->boundary)--;
      t->vertices[j] = v;
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
        e->incident_triangles[1] = t;
      }
      else {
        e = &(m->edges[edge_count]);
        edge_count++;
        e->incident_triangles[0] = t;
        e->incident_triangles[1] = NULL;
        add_incident_vertices_and_edge(v1, v2, e);
      }
      opposite_edge = (j+2) % 3;
      t->edges[opposite_edge] = e;
    }
  }
  return edge_count;
}


/* Finds all boundary edges, and creates the cycles of
 * edges for each boundary component.
 */
void calc_boundaries(mesh *m)
{
  int i;
  std::list<edge *>* cycle;
  edge* e;
  edge* e_start;
  int* track_boundary_edges = (int*) malloc(m->ranks[1] * sizeof(int));
  m->boundary_edges = (int*) malloc(m->ranks[1] * sizeof(int));
  if (!track_boundary_edges || !m->boundary_edges) {
    printf("Memory allocation error in calc_boundaries.\n");
    exit(1);
  }
  m->boundary_count=0;
  for (i=0; i<m->ranks[1]; i++) {
    if (m->edges[i].incident_triangles[1] == NULL) {
      m->boundary_edges[i] = 1;
      track_boundary_edges[i] = 1;
    } else {
      m->boundary_edges[i] = 0;
      track_boundary_edges[i] = 0;
    }
  }
  while((i=get_next_component_start(track_boundary_edges, m->ranks[1])) >= 0) {
    track_boundary_edges[i] = 0;
    cycle = &(m->boundary_cycles[m->boundary_count]);
    (m->boundary_count)++;
    e_start = &(m->edges[i]);
    (*cycle).push_front(e_start);
    e = e_start;
    while ((e=get_next_boundary_edge(e)) != e_start) {
      (*cycle).push_back(e);
      track_boundary_edges[e->index] = 0;
    }
  }
  free(track_boundary_edges);
}

/* Returns the index of the first non-zero entry in
 * track_boundary_edges. This corresponds with the first
 * boundary edge not yet put in a boundary cycle.
 * Returns -1 if all entries are 0.
 */
int get_next_component_start(int* track_boundary_edges, int n)
{
  int i;
  for (i=0; i<n; i++) {
    if (track_boundary_edges[i] != 0) {
      return i;
    }
  }
  return -1;
}

/* Assuming that edge e is a boundary edge, 
 * function returns the next edge if following along
 * the induced orientation along the boundary.
 */
edge* get_next_boundary_edge(edge* e)
{
  vertex* v = (vertex*)e->vertices[1];
  return (edge*)v->incident_edges[0];
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
  vertex* v;
  vertex* iv[max_degree * 2 + 1];
  edge* ie[max_degree * 2 + 1];
  int head, tail, i;
  
  for (i=0; i < m->ranks[0]; i++) {
    v = &(m->vertices[i]);
    head = max_degree;
    tail = max_degree;
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

/* Constructs the list of link edges around each vertex, with ordering agreeing
 * with the existing cyclic ordering of the incident vertices, edges and
 * triangles.
 *
 * Note that a link edge is not incident to the vertex, but shares a common
 * triangle (so the union of all link edges of a vertex is a circle around the
 * vertex).
 */
void add_link_edges(mesh *m)
{
  int i, j, k;
  vertex *v;
  triangle *t;
  for (i=0; i < m->ranks[0]; i++) {
    v = &(m->vertices[i]);
    for (j=0; j < v->degree - v->boundary; j++) {
      t = (triangle *)(v->incident_triangles[j]);
      for (k=0; k<3; k++) {
        if (t->vertices[k] == v) {
          v->link_edges[j] = t->edges[k];
        }
      }
    }
  }
}

void clear_geodesic_lists(mesh* m)
{
  vertex* v;
  geodesic g(0,0,0,0);

  for (int i=0; i<m->ranks[0]; ++i) {
    v = &(m->vertices[i]);
    v->geodesics.clear();
  }
}

/* Links the hessian data at each vertex to the hessian data recorded
 * at each triangle (the latter is easier to calculate, the former
 * is easier to access when performing Ricci flow).
 */
void construct_vertex_hessian_pointers(mesh *m)
{
  vertex *v;
  triangle *t;
  int i, j, k, vi;
  for (vi=0; vi < m->ranks[0]; vi++) {
    v = &(m->vertices[vi]);
    for (i=0; i < v->degree - v->boundary; i++) {
      t = (triangle *)(v->incident_triangles[i]);
      k = get_vertex_position_in_triangle(v,t);
      for (j=0; j<3; j++) {
        v->dtheta_du[i][j] = &(t->hessian[j][k]);
      }
    }
  }
}

/* constructs an array of pointers for each
 * vertex to the inner angles at the vertex.
 * The ordering of the array matches the ordering
 * of the incident triangles.
 */
void construct_vertex_inner_angle_pointers(mesh *m)
{
  int i, j, k;
  vertex *v;
  triangle *t;
  for(i=0; i<m->ranks[0]; i++) {
    v = &(m->vertices[i]);
    for (j=0; j < v->degree-v->boundary; j++) {
      t = (triangle*)v->incident_triangles[j];
      k = get_vertex_position_in_triangle(v, t);
      v->inner_angles[j] = &(t->inner_angles[k]);
    }
  }
}

/* Assuming that vertex v is incident to triangle t, function
 * returns the order (0, 1 or 2) of the vertex within the triangle.
 *
 * Exits rudely if the vertex is not incident.
 */
int get_vertex_position_in_triangle(vertex *v, triangle *t)
{
  int i;
  for (i=0; i<3; i++) {
    if (v == (t->vertices[i])) {
      return i;
    }
  }
  printf("Vertex is not incident to triangle - something is wrong!\n");
  exit(1);
}


/* Assuming that vertex v is incident to edge e,
 * functin returns a pointer to the other vertex
 * incident to e.
 *
 * Exits rudely if the vertex is not incident.
 */
vertex* get_other_vertex(edge* e, vertex* v)
{
  if (e->vertices[0] == v) return e->vertices[1];
  if (e->vertices[1] == v) return e->vertices[0];
  printf("Vertex is not incident to edge - error in get_other_vertex.\n");
  exit(1);
}
/* Assuming that edge e is incident to vertex v, function
 * returns the index of the edge within v->incident_edges.
 *
 * Exits rudely if the edge is not incident.
 */
int get_edge_position_at_vertex(edge* e, vertex* v)
{
  for (int i=0; i<v->degree; ++i) {
    if (e == (v->incident_edges[i])) {
      return i;
    }
  }
  printf("Edge is not incident to vertex - something is wrong!\n");
  exit(1);
}

bool edges_share_triangle(edge* e1, edge* e2)
{
  for (int i=0; i<2; ++i) {
    triangle* t = e1->incident_triangles[i];
    if (!t) return false;
    for (int j=0; j<3; ++j) {
      if (t->edges[j] == e2) return true;
    }
  }
  return false;
}

/* Finds the edge incident to v, immediately clockwise from ie[0],
 * and records it in ie[-1]. Does the same for the vertices.
 * If there is no edge clockwise from ie[0] (it is on the boundary)
 * then the function returns 0. Otherwise, it returns 1.
 */
int find_clockwise_edge(vertex* v, edge* ie[], vertex* iv[])
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
int find_counterclockwise_edge(vertex* v, edge* ie[], vertex* iv[])
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
void sort_incident_edges_and_vertices(vertex* v, edge* ie[], vertex* iv[])
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
triangle* get_common_triangle(edge *e1, edge *e2)
{
  int i, j;
  triangle *t;
  for (i=0; i<2; i++) {
    t = (triangle *)(e1->incident_triangles[i]);
    if (t) {
      for (j=0; j<3; j++) {
        if (t->edges[j] == e2)
          return t;
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

/* Changes the s parameters at all vertices to the passed values.
 * Also calculates f at the new parameter values.
 * Since f is defined by an integral which can only
 * be numerically estimated, the updates to s occurs
 * in stages to allow for Simpson's rule to be used.
 *
 * n is the goal number of divisions for the integration
 * domain.
 *
 * Function returns the updated value for f(s).
 */
double update_f_and_s(mesh *m,  double *s, int n)
{
  int i;
  double s1, s0;
  double ds;
  vertex *v;
  for (i=0; i < m->ranks[0]; i++) {
    v = &(m->vertices[i]);
    s0 = v->s;
    s1 = s[i];
    ds = fabs(s1 - s0) / n;
    m->f += integrate(&curvature_integrand, s0, s1, ds, (void *)v);
  }
  return m->f;
}

/* Updates the radius parameters of all vertices to the 
 * passed values. Also recalculates the edge lengths
 * using the new radii.
 */
void update_s_and_edge_lengths(mesh *m, MV_Vector_double &s)
{
  int i;
  for (i=0; i < m->ranks[0]; i++) {
    m->vertices[i].s = s[i];
  }
  calc_edge_lengths(m);
}

/* Callback function used to integrate curvatures. After calling,
 * the s parameter of the vertex will be updated to s, and the lengths
 * of all incident edges will be recalculated accordingly.
 *
 * Integrand is K / s because of change of variables u -> s = e^u.
 */
double curvature_integrand(double s, void *instance)
{
  int i;
  edge *e;
  vertex *v = (vertex *)instance; 
  v->s = s;
  for (i=0; i < v->degree; i++) {
    e = (edge *)(v->incident_edges[i]);
    calc_edge_length(e);
  }
  return calc_curvature(v) / s;
}

/* Calculates the curvatures at all vertices and stores
 * the values in the array K.
 *
 * TODO: Update inner angles when edge lengths are calculated.
 */
void calc_curvatures(mesh *m, MV_Vector_double &K)
{
  int i;
  vertex *v;
  calc_inner_angles(m);
  for (i=0; i < m->ranks[0]; i++) {
    v = &(m->vertices[i]);
    K[i] = calc_curvature(v);
  }
}

/* Calculates interior angles in the triangles
 * of the mesh (using hyperbolic trig).
 */
void calc_inner_angles(mesh *m)
{
  triangle *t;
  edge *e;
  int i, j;
  double cl[3];
  double sq[3];
  double cos_angle;
  for(i=0; i < m->ranks[2]; i++) {
    t = &(m->triangles[i]);
    for (j=0; j<3; j++) {
      cl[j] = ((edge *)(t->edges[j]))->cosh_length_minus1;  
      sq[j] = cl[j] * cl[j];
    }
    for (j=0; j<3; j++) {
      cos_angle = cl[(j+2)%3] + cl[(j+1)%3] + cl[(j+2)%3] * cl[(j+1)%3] - cl[j];
      cos_angle = cos_angle / sqrt((2 * cl[(j+2)%3] + sq[(j+2)%3]) 
          * (2 * cl[(j+1)%3] + sq[(j+1)%3]));
      if (cos_angle > 1.0) {
        cos_angle = 1.0;  // Protect against rounding errors.
      }
      t->inner_angles[j] = acos(cos_angle);
      t->cos_angles[j] = cos_angle;
      t->sin_angles[j] = sin(t->inner_angles[j]);
#ifdef DEBUG_MODE
      if (!is_number(cos_angle) || !is_number(t->inner_angles[j]) || !is_number(t->sin_angles[j])) {
        printf("NaN found in calc_inner_angles.\n");
      }
#endif
    }
  }
}

/* Calculates the curvature K at vertex v,
 * and returns the result. Note that
 *      K = 2Pi - (1+b)Sum(interior angles )
 * where the sum is over all triangles incident to v, with
 * b = 1 if v is a boundary vertex, 0 otherwise.
 */
double calc_curvature(vertex *v)
{
  int i, j;
  double angle_sum = 0;
  triangle *t;
  for (i=0; i < v->degree - v->boundary; i++) {
    t = (triangle *)(v->incident_triangles[i]);
    j = get_vertex_position_in_triangle(v, t);
    angle_sum += t->inner_angles[j];
  }
  return (2 - v->boundary) * M_PI - angle_sum;
}

void calc_edge_lengths(mesh *m) 
{
  int i;
  double min_length = 100;
  for (i=0; i < m->ranks[1]; i++) {
    double length = calc_edge_length(&(m->edges[i]));
    if (min_length > length) {
      min_length = length;
    }
#ifdef DEBUG_MODE
    if (length < edge_length_threshold) {
      printf("Short edge detected.\n");
    }
#endif
  }
  if (min_length < edge_length_threshold) {
    // throw min_length;
  }
}
double calc_edge_length (edge *e)
{
  double a, b, c, aa, bb;
  a = ((vertex *)(e->vertices[0]))->s;
  b = ((vertex *)(e->vertices[1]))->s;
  c = e->cos_angle;
  aa = a * a;
  bb = b * b;
  e->cosh_length_minus1 = (2 * (aa + bb) - 4*a*b*c) / 
    ((1 - aa)*(1 - bb));
  e->sinh_length = sqrt(2 * e->cosh_length_minus1 
      + e->cosh_length_minus1 * e->cosh_length_minus1);
#ifdef DEBUG_MODE
  if ((!is_number(e->cosh_length_minus1)) || (!is_number(e->sinh_length))) {
    printf("NaN found in calc_edge_length.\n");
  }
#endif
  return e->cosh_length_minus1;
}

void calc_boundary_lengths(mesh* m)
{
  int i;
  double length;
  std::list<edge*>::iterator j;

  for (i=0; i<m->boundary_count; i++) {
    length = 0;
    j = m->boundary_cycles[i].begin();
    for(; j != m->boundary_cycles[i].end(); j++) {
      length += acosh(1.0 + (*j)->cosh_length_minus1);
    }
    m->boundary_lengths[i] = length;
  }
}

/* Look for edges shorter than edge_length_threshold and
 * collapse them.
 */
void collapse_short_edges(mesh* m)
{
  int* vertex_map = new int[m->ranks[0]];
  int** inverse_map;
  double* s = new double[m->ranks[0]];
  int** edge_incidences;
  double** edge_angles;
  int original_number_of_vertices = m->ranks[0];
  initialize_double_arrays_of_edges(&edge_incidences, &edge_angles,
      &inverse_map, original_number_of_vertices);
  record_edge_angles(m, edge_incidences, edge_angles);
  record_vertex_radii(m, s);
  construct_vertex_equivalence_relation(m, vertex_map);
  construct_inverse_vertex_map(m, vertex_map, inverse_map);
  remap_collapsed_triangles(m, vertex_map); 

  construct_simplices(m, m->fd);
  sort_cyclic_order_at_vertices(m);
  add_wings_to_bivalent_vertices(m, edge_incidences, edge_angles);
  remove_isolated_vertices(m, m->fd);
  sort_cyclic_order_at_vertices(m);
  add_link_edges(m);
  clear_boundary_data(m);
  calc_boundaries(m);
  construct_vertex_hessian_pointers(m);
  construct_vertex_inner_angle_pointers(m);
  copy_radii_and_angle_data(m, inverse_map, s, edge_incidences, edge_angles);
  delete [] vertex_map;
  delete [] s;
  deallocate_double_arrays_of_edges(edge_incidences, edge_angles, 
      inverse_map, original_number_of_vertices);
}

void add_wings_to_bivalent_vertices(mesh* m, int** edge_incidences, 
    double** edge_angles)
{
  face* face_node = m->fd->face_head;
  while (face_node->next != NULL) {
    face_node = face_node->next;
  }
  for (int i=0; i<m->ranks[0]; ++i) {
    vertex* v = &(m->vertices[i]);
    if (v->degree == 2) {
      int vi[4];
      get_boundary_adjacent_vertices(v, vi);
      ++(m->vertices[vi[0]].degree); // to prevent nearby bivalent conflicts
      ++(m->vertices[vi[3]].degree);
      face_node = add_face_node(face_node, v->index, vi[0], vi[1]);
      face_node = add_face_node(face_node, v->index, vi[2], vi[3]);
      m->fd->number_of_triangles += 2;
      add_pi_edge(edge_incidences, edge_angles, v->index, vi[0]);
      add_pi_edge(edge_incidences, edge_angles, v->index, vi[3]);
    }
  }
  construct_simplices(m, m->fd);
  sort_cyclic_order_at_vertices(m);
}

void get_boundary_adjacent_vertices(vertex* v, int* vi)
{
  vertex* v1 = v->incident_vertices[0];
  vertex* v2 = v->incident_vertices[1];
  vi[1] = v1->index;
  vi[2] = v2->index;
  vi[0] = v1->incident_vertices[0]->index;
  vi[3] = v2->incident_vertices[v2->degree - 1]->index;
}

void add_pi_edge(int** ei, double** ea, int v0, int v1)
{
  int i=0;
  while (ei[v0][i] != -1) ++i;
  ei[v0][i] = v1;
  ea[v0][i] = -1;
  i=0;
  while (ei[v1][i] != -1) ++i;
  ei[v1][i] = v0;
  ea[v1][i] = -1;
}

void clear_boundary_data(mesh* m)
{
  free(m->boundary_edges);
  for (int i=0; i<m->boundary_count; ++i) {
    m->boundary_cycles[i].clear();
  }
}
void initialize_double_arrays_of_edges(int*** ei, double*** ea, 
    int*** im, int n)
{
  *ei = new int*[n];
  *ea = new double*[n];
  *im = new int*[n];
  for (int i=0; i<n; ++i) {
    (*ei)[i] = new int[max_degree];
    for (int j=0; j<max_degree; ++j) {
      (*ei)[i][j] = -1;
    }
    (*ea)[i] = new double[max_degree];
    (*im)[i] = new int[max_degree];
  }
}

void deallocate_double_arrays_of_edges(int** ei, double** ea, 
    int** im, int n)
{
  for (int i=0; i<n; ++i) {
    delete [] ei[i];
    delete [] ea[i];
    delete [] im[i];
  }
  delete [] ei;
  delete [] ea;
  delete [] im;
}

void record_edge_angles(mesh* m, int** ei, double** ea)
{
  for(int i=0; i<m->ranks[0]; ++i) {
    vertex* v = &(m->vertices[i]);
    for(int j=0; j<v->degree; ++j) {
      vertex* v1 = v->incident_vertices[j];
      edge* e = v->incident_edges[j];
      ei[i][j] = v1->index;
      ea[i][j] = e->cos_angle;
    }
    ei[i][v->degree] = -1;
  }
}

void record_vertex_radii(mesh* m, double* s)
{
  for(int i=0; i<m->ranks[0]; ++i) {
    vertex* v = &(m->vertices[i]);
    s[i] = v->s;
  }
}

void construct_vertex_equivalence_relation(mesh* m, int* vertex_map)
{
  for (int i=0; i<m->ranks[0]; ++i) {
    vertex_map[i] = -1;
  }
  for (int i=0; i<m->ranks[0]; ++i) {
    vertex* v = &(m->vertices[i]);
    if (vertex_map[i] == -1) {
      vertex_map[i] = i;
    }
    for (int j=0; j<v->degree; ++j) {
      edge* e = v->incident_edges[j];
      if (e->cosh_length_minus1 < edge_length_threshold) {
        vertex* opposite_v = get_other_vertex(e, v);
        if (opposite_v->index > v->index) {
          vertex_map[opposite_v->index] = vertex_map[v->index];
        }
      }     
    }
  }
}

void construct_inverse_vertex_map(mesh* m, int* vertex_map, int** inverse_map)
{
  int position = 0;
  for (int i=0; i<m->ranks[0]; ++i) {
    inverse_map[i][0] = 0;  // this slot records the number of inverse points
  }
  for (int i=0; i<m->ranks[0]; ++i) {
    if (vertex_map[i] == i) {
      int count = 0;
      for (int j = i; j<m->ranks[0]; ++j) {
        if (vertex_map[j] == i) {
          inverse_map[position][count+1] = j;
          ++count;
        }
      }
      inverse_map[position][0] = count;
      ++position;
    }
  }
}

void remap_collapsed_triangles(mesh* m, int* vertex_map)
{
  face* face_head = (face*)malloc(sizeof(face));
  face_head->next = NULL;
  face* face_node = face_head;
  int number_of_triangles = 0;
  for (int i=0; i<m->ranks[2]; ++i) {
    triangle* t = &(m->triangles[i]);
    int count = count_collapsed_edges(t);
#ifdef DEBUG_MODE
    if (count == 2) {
      printf("Only two collapsed edges in a triangle!\n");
    }
#endif
    if (count == 0) {
      ++number_of_triangles;
      int v[3];
      for (int j=0; j<3; ++j) {
        v[j] = vertex_map[t->vertices[j]->index];
      }
      face_node = add_face_node(face_node, v[0], v[1], v[2]); 
    }
  }
  deallocate_face(m->fd->face_head);
  m->fd->face_head = face_head;
  m->fd->number_of_triangles = number_of_triangles;
}

int count_collapsed_edges(triangle* t)
{
  int count = 0;
  for (int i=0; i<3; ++i) {
    edge* e = t->edges[i];
    if (e->cosh_length_minus1 < edge_length_threshold) {
      ++count;
    }
  }
  return count;
}

void copy_radii_and_angle_data(mesh* m, int** inverse_map, 
    double* s, int** edge_incidences, double** edge_angles)
{
  for (int i=0; i<m->ranks[0]; ++i) {
    vertex* v = &(m->vertices[i]);
    int index = inverse_map[i][1];
    if (s[index] > radii_threshold) {
      v->s = s[index];
    } else {
      v->s = radii_threshold;
    }
  }
  for (int i=0; i<m->ranks[1]; ++i) {
    edge* e = &(m->edges[i]);
    int i0 = e->vertices[0]->index;
    int i1 = e->vertices[1]->index;
    e->cos_angle = find_inverse_edge_angle(i0, i1, inverse_map,
        edge_incidences, edge_angles);
  }
}

double find_inverse_edge_angle(int i0, int i1, int** inverse_map,
    int** ei, double** ea)
{
  double angle = 0;
  double min_angle = 1;
  int n = inverse_map[i0][0];
  int m = inverse_map[i1][0];
  int* list0 = &(inverse_map[i0][1]);
  int* list1 = &(inverse_map[i1][1]);
  for (int i=0; i<n; ++i) {
    for (int j=0; j<m; ++j) {
      int k0 = list0[i];
      int k1 = list1[j];
      int pos = find_old_edge_position(ei, k0, k1);
      if (pos != -1) {
        angle = ea[k0][pos];
        if (angle < min_angle) {
          min_angle = angle;
        }
      }
    }
  }
  if (min_angle > 0) {
    printf("Error finding edge angle.\n");
  }
  return min_angle;
}

int find_old_edge_position(int** ei, int k0, int k1)
{
  for (int i=0; i<max_degree; ++i) {
    if (ei[k0][i] == k1) return i;
    if (ei[k0][i] == -1) return -1;
  }
}
/* Calculates the the negative of the Hessian matrix [dK_i / du_j], given 
 * that edge angles and radius parameters are set.
 */
void calc_hessian(mesh *m)
{
  double dtheta_dl[3][3];
  double dl_ds[3][3];
  double ds_du[3][3];
  double entry, entry0, entry1;
  triangle *t;
  int i, j, k, l, k0, k1, count;
  int nz;
  vertex *v0, *v1;
  edge *e;

  for (i=0; i < m->ranks[2]; i++) {
    t = &(m->triangles[i]);
    calc_dtheta_dl(t, dtheta_dl);
    calc_dl_ds(t, dl_ds);
    calc_ds_du(t, ds_du);
    calc_3_matrix_product(dtheta_dl, dl_ds, ds_du,
        t->hessian);
  }


  nz = m->ranks[0] + 2 * m->ranks[1];
  double *val = (double *)malloc(nz * sizeof(double));
  int *row_ind = (int *)malloc(nz * sizeof(int));
  int *col_ind = (int *)malloc(nz * sizeof(int));
  if (!val || !row_ind || !col_ind) {
    printf("Memory allocation error in calc_hessian.\n");
    exit(1);
  }
  count = 0;
  for (i=0; i < m->ranks[0]; i++) {
    v0 = &(m->vertices[i]);
    entry = 0;
    for (j=0; j < v0->degree - v0->boundary; j++) {
      t = (triangle *)(v0->incident_triangles[j]);
      k = get_vertex_position_in_triangle(v0, t);
      entry += *(v0->dtheta_du[j][k]);
    }
    // entry *= 1 + v0->boundary; // boundary vertices count double
    val[count] = entry;
    row_ind[count] = i;
    col_ind[count] = i;
    count++;
  }
  for (i=0; i < m->ranks[1]; i++) {
    e = &(m->edges[i]);
    v0 = (vertex *)(e->vertices[0]);
    v1 = (vertex *)(e->vertices[1]);
    entry0 = 0; entry1 = 0;
    for (j=0; j<2 && e->incident_triangles[j]; j++) {
      t = (triangle *)(e->incident_triangles[j]);
      k0 = get_vertex_position_in_triangle(v0, t);
      k1 = get_vertex_position_in_triangle(v1, t);
      entry0 += t->hessian[k0][k1];
      entry1 += t->hessian[k1][k0];
    }
    val[count] = entry0;
    row_ind[count] = v0->index;
    col_ind[count] = v1->index;
    count++;
    val[count] = entry1;
    row_ind[count] = v1->index;
    col_ind[count] = v0->index;
    count++;
  }
  Coord_Mat_double A(m->ranks[0], m->ranks[0], nz, val, row_ind, col_ind);
  m->hessian = A;

  free(val);
  free(row_ind);
  free(col_ind);

}

int calc_number_of_nonzero_hessian_entries(mesh *m)
{
  int i;
  vertex v;
  int count = 0;
  for (i=0; i < m->ranks[0]; i++) {
    v = (m->vertices[i]);
    count += 3 * (v.degree - v.boundary);
  }
  return count;
}

/* Calculate the matrix of derivatives [d theta/dl],
 * and stores the result in A.
 */
void calc_dtheta_dl(triangle *t, double A[3][3])
{
  double d1, d2;
  double cosh_l[3], sinh_l[3];
  int i, j, k, j1, j2;
  for (i=0; i<3; i++) {
    cosh_l[i] = 1.0 + ((edge *)t->edges[i])->cosh_length_minus1;
    sinh_l[i] = ((edge *)t->edges[i])->sinh_length;
  }
  for (i=0; i<3; i++) {
    A[i][i] = sinh_l[i] / (t->sin_angles[i] *
        sinh_l[(i+1)%3] * 
        sinh_l[(i+2)%3]);
    for (j1=1; j1<3; j1++) {
      j2 = 3 - j1;
      j = (i+j1) % 3;
      k = (i+j2) % 3; // i, j, k all distinct
      d1 = cosh_l[j] * t->cos_angles[i] / 
        sinh_l[j];
      d2 = cosh_l[k] / sinh_l[k];
      A[i][j] = (d1 - d2) / t->sin_angles[i];
#ifdef DEBUG_MODE
      if (!is_number(d1) || !is_number(d2) || !is_number(A[i][j]) || !is_number(A[i][i])) {
        printf("NaN found in calc_dtheta_dl.\n");
      }
#endif
    }
  }
}

/* Calculate the matrix of derivatives [dl/ds],
 * and stores the result in A.
 */
void calc_dl_ds(triangle *t, double A[3][3])
{
  int i, j, k, j1, j2;
  double d1, d2;
  double s[3], cos_e[3], cosh_l[3], sinh_l[3];
  for (i=0; i<3; i++) {
    s[i] = ((vertex *)t->vertices[i])->s;
    cosh_l[i] = 1.0 + ((edge *)t->edges[i])->cosh_length_minus1;
    sinh_l[i] = ((edge *)t->edges[i])->sinh_length;
    cos_e[i] = ((edge *)t->edges[i])->cos_angle;
  }
  for (i=0; i<3; i++) {
    A[i][i] = 0;
    for (j1=1; j1<3; j1++) {
      j2 = 3 - j1;
      j = (i+j1) % 3;
      k = (i+j2) % 3; // i, j, k all distinct
      d1 = (2*s[j] * (1 + s[k]*s[k]) - 4*s[k] * cos_e[i] ) /
             ((1 - s[j]*s[j]) * (1 - s[k]*s[k]));
      d2 = 2 * s[j] * cosh_l[i] / (1 - s[j]*s[j]);
      A[i][j] = (d1 + d2) / sinh_l[i]; 
#ifdef DEBUG_MODE
      if (!is_number(d1) || !is_number(d2) || !is_number(A[i][j])) {
        printf("NaN found in calc_dl_ds.\n");
      }
#endif
    }
  }
}

/* Calculate the matrix of derivatives [ds/du],
 * and stores the result in A.
 */
void calc_ds_du(triangle *t, double A[3][3])
{
  int i, j;
  for (i=0; i<3; i++) {
    A[i][i] = ((vertex *)t->vertices[i])->s;
    for (j=1; j<3; j++) {
      A[i][(i+j)%3] = 0;
    }
  }
}

/* Calculate the product ABC, and stores the result in D.
 */
void calc_3_matrix_product(double A[3][3], double B[3][3],
    double C[3][3], double D[3][3])
{
  int i, j, k, l;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      D[i][j] = 0;
      for (k=0; k<3; k++) {
        for (l=0; l<3; l++) {
          D[i][j] += A[i][k] * B[k][l] *
            C[l][j];
        }
      }
    }
  }
}

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
  v1->incident_vertices[v1->degree] = v2;
  v1->incident_edges[v1->degree] = e;
  (v1->degree)++;;
  (v1->boundary)++;
  v2->incident_vertices[v2->degree] = v1;
  v2->incident_edges[v2->degree] = e;
  (v2->degree)++;;
  (v2->boundary)++;
  e->vertices[0] = v1;
  e->vertices[1] = v2;
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

double min(double *x, int n)
{
  int i;
  double smallest = x[0];
  for (i=0; i<n; i++) {
    smallest = x[i] < smallest ? x[i] : smallest;
  }
  return smallest;
}

double max(double *x, int n)
{
  int i;
  double largest = x[0];
  for (i=0; i<n; i++) {
    largest = x[i] > largest ? x[i] : largest;
  }
  return largest;
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
  printf("]  LE:[");
  for (j=0; j < v->degree - v->boundary; j++) {
    e = (edge *)(v->link_edges[j]);
    k= e->index;
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
    printf("V:[%i, %i]  T:[%i]", v1, v2, t1);
  }
  else {
    t2 = ((triangle *)t)->index;
    printf("V:[%i, %i]  T:[%i, %i]", v1, v2, t1, t2);
  }
  printf(" (Cos(A), Cosh(L)) = (%f, %f)\n", e->cos_angle, 1.0 + e->cosh_length_minus1);
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

/* Function should return true is x is a (possibly infinite)
 * floating point number, return false if it a nan.
 */
bool is_number(double x)
{
  return (x == x);
}
