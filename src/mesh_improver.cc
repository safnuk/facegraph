
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <list>

#include <lbfgs.h>
#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"
#include "cfile_io.h"
#include "geodesic.h"
#include "cmesh.h"
#include "ccirclepack.h"

#include "mesh_improver.h"

using namespace std;

/* Code looks for vertices which prevent Ricci flow convergence
 * and subdivides the adjacent triangles to increase the triangle
 * edge angles.
 *
 * Algorithm reconstructs the mesh and finds the new circle packing metric.
 */
int find_problem_vertices(mesh* m, filedata* fd)
{
  int count = 0;
  int trouble_count = 0;
  int new_middle_points[3];
  int incident_points[3];
  int trivalent_point;
  int* triangles_to_keep = new int[m->ranks[2]];
  for (int i=0; i<m->ranks[2]; ++i) {
    triangles_to_keep[i] = 1;
  }
  int * vertex_hit_list = new int[m->ranks[0]];
  for (int i=0; i<m->ranks[0]; ++i) {
    vertex_hit_list[i] = 0;
  }
  face* face_head = (face*)malloc(sizeof(face));
  face_head->next = NULL;
  face* face_node = face_head;
  for (int i=0; i<m->ranks[0]; ++i) {
    double angle_sum = 0;
    vertex* v = &(m->vertices[i]);
    if (v->degree + v->boundary <= 4) {
      ++count;
      for (int j=0; j<v->degree - v->boundary; ++j) {
        angle_sum += acos(v->link_edges[j]->cos_angle);
      }
      angle_sum *= 1 + v->boundary;
      if ((angle_sum <= 2*M_PI + 0.2) && (!vertex_hit_list[v->index])) {
        ++trouble_count;
        list_incident_vertices(v, incident_points, vertex_hit_list);
        trivalent_point = v->index;
        remove_triangles_from_keep_list(v, triangles_to_keep);
        construct_incident_edge_bisectors(m, fd, v, new_middle_points);
        face_node = add_new_triangles_for_trivalent(face_node, trivalent_point,
            incident_points, new_middle_points, v->degree, v->boundary);
      }
    }
  }
  construct_new_triangle_list(fd, face_head, triangles_to_keep);
  deallocate_mesh(m);
  initialize_mesh(m, fd);
  printf("Improved %i vertices (out of %i possible problems).\n",
    trouble_count, count);
  if (trouble_count > 0) {
    calc_circlepack_metric(m);
  }
  delete [] triangles_to_keep;
  delete [] vertex_hit_list;
  return trouble_count;
}

/* Function lists the indices of the vertices incident to v,
 * and stores the list in incident_points.
 * Also records these vertices in vertex_hit_list, as after
 * the subdivision they will have an additional incident edge
 * (and hence no longer be a problem vertex if they were before).
 */
void list_incident_vertices(vertex* v, int* incident_points, 
    int* vertex_hit_list)
{
  for (int i=0; i<v->degree; ++i) {
    vertex* v2 = v->incident_vertices[i];
    incident_points[i] = v2->index;
    vertex_hit_list[v2->index] = 1;
  }
}

/* Records the triangles incident to v, as they will be subdivided by the
 * algorithm, hence not longer needed.
 */
void remove_triangles_from_keep_list(vertex* v, int* triangles_to_keep)
{
  for (int i=0; i<v->degree - v->boundary; ++i) {
    triangle* t = v->incident_triangles[i];
    triangles_to_keep[t->index] = 0;
  }
}

/* Finds the midpoints of all edges incident to v, and adds the
 * resulting coordinates in the point list of fd. The indices of the
 * new vertices are recorded in new_middle_points.
 */
void construct_incident_edge_bisectors(mesh* m, filedata* fd, vertex* v, 
                                       int* new_middle_points)
{
  _point* point_node = fd->point_head;
  while (point_node->next != NULL){
    point_node = point_node->next;
  }
  for (int i=0; i<v->degree; ++i) {
    edge* e = v->incident_edges[i];
    int j0 = e->vertices[0]->index;
    int j1 = e->vertices[1]->index;
    double x = (m->coordinates[j0].x + m->coordinates[j1].x) / 2;
    double y = (m->coordinates[j0].y + m->coordinates[j1].y) / 2;
    double z = (m->coordinates[j0].z + m->coordinates[j1].z) / 2;
    new_middle_points[i] = fd->number_of_points;
    point_node = add_point_node(point_node, x, y, z);
    ++(fd->number_of_points);
  }
}

/* Lists out the new triangles formed from subdivision (each
 * original incident triangle is divided into three triangles).
 */
face* add_new_triangles_for_trivalent(face* face_node, int p, int* v, int* m, 
    int degree, int boundary)
{
 for (int i=0; i<degree-boundary; ++i) {
    face_node = add_face_node(face_node, p, m[i], m[(i+1)%degree]);
    face_node = add_face_node(face_node, m[i], v[i], v[(i+1)%degree]);
    face_node = add_face_node(face_node, m[i], v[(i+1)%degree], 
        m[(i+1)%degree]);
 } 
 return face_node;
}
