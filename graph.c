// graph.c
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>
#include <iostream>

#include "coord_double.h"
#include "comprow_double.h"
#include "compcol_double.h"
#include "mvblasd.h"

#include "cfile_io.h"
#include "geodesic.h"
#include "cmesh.h"
#include "graph.h"

/* Finds the boundary nearest to each vertex and
 * calculates the geodesic realizing the shortest
 * distance.
 * 
 * Algorithm does not properly handle graphs with
 * edges which have the same boundary on both sides.
 */
void calc_vertex_boundary_distances(mesh *m)
{
        std::list<vertex*> active;
        int i;
        for (i=0; i<m->boundary_count; i++) {
                clear_geodesic_lists(m);
                create_active_list(m, i, active);
                run_through_active_list(m, i, active);
        }
}

/* Initializes the list of currently active vertices
 * with the vertices along boundary cycle b.
 */
void create_active_list(mesh *m, int b, std::list<vertex*>& active)
{
        vertex* v;
        edge* e;
        std::list<edge*>::iterator cycle; 
        double position = 0;
        for (cycle=m->boundary_cycles[b].begin(); cycle!=m->boundary_cycles[b].end(); cycle++) {
                e = *cycle;
                v = (vertex*)(e->vertices[0]);
                v->shortest_paths[b].assign_values(b, 0, position, M_PI_2);
                calc_vertex_config(v, v->shortest_paths[b], true);
                active.push_back(v);
                position += acosh(e->cosh_length);
        }
        // TODO: remove
        int length = active.size();
        return;
}

/* For every vertex in the active list, calculates the geodesic
 * distances to boundary b for every vertex incident to it and
 * further away from the boundary. If such a vertex has been hit
 * with all incident vertices which are nearer to the boundary,
 * calculate the average of all these hits and add  the vertex
 * to the active list.
 *
 * Continue until no more vertices are active.
 */
void run_through_active_list(mesh *m, int b, std::list<vertex*>& active)
{
        std::list<vertex*>::iterator i = active.begin();
        vertex* v;
        vertex* v1;
        geodesic g;
        while (i!=active.end()) {
                v = *i;
                std::list<int>::iterator j = v->vc.farther_vertices.begin();
                for (;j!=v->vc.farther_vertices.end(); j++) {
                        calc_next_vertex_geodesic(m, v, *j, g, b, NULL);
                        v1 = (vertex*)(v->incident_vertices[*j]);
                        add_geodesic_to_vertex(v1, v, g);
                        if((v1->vc.closer_vertices.empty()) && (v1->shortest_paths[b].boundary == -1)) {
                                geodesic error = calc_average_geodesic(v1);
                                if (error > kErrorThreshold) { // Non-optimal geodesics in the list
                                        if (recalc_vertex_geodesic(v1) == kVertexActive) {
                                                active.push_back(v1);
                                        }
                                } else {
                                        active.push_back(v1);
                                }
                        }
                }
                i = active.erase(i);
        }
}

/* For a vertex which has been hit with non-optimal geodesics,
 * this function recalculate the shortest path by restricting 
 * to optimal geodesics.
 *
 * Returns the new status of the vertex (active or not).
 */
int recalc_vertex_geodesic(vertex* v) {
        geodesic g_min = *(std::min_element(v->geodesics.begin(), 
                                v->geodesics.end()));
        calc_vertex_config(v, g_min, false);
        std::list<geodesic>::iterator i = v->geodesics.begin();
        while (i != v->geodesics.end()) {
                std::list<*vertex>::iterator j = find(v->vc.closer_vertices.begin(),
                                v->vc.closer_vertices.end(), (*i).originating_vertex);
                if (j == v->vc.closer_vertices.end()) {  // vertex not in closer list
                        i = v->geodesics.remove(i);
                } else {
                        v->vc.closer_vertices.remove(j);
                        ++i;
                }
        }
        if (v->vc.closer_vertices.empty()) {
                geodesic error = calc_average_geodesic(v);
                if (error > kErrorThreshold) {
                        printf("We're in trouble!\n");
                }
                return kVertexActive;
        } else {
                return kVertexNotActive;
        }
}

/* Calculates which incident vertices are further away from v 
 * and which are closer. The calculation is simpler for
 * points known to be on the boundary, hence the flag on_boundary. 
 */
void calc_vertex_config(vertex* v, geodesic const& g, bool on_boundary)
{
        const double fudge_factor = 1.0000000001;
        int i;
        geodesic g_next;
        v->vc.farther_vertices.clear();
        v->vc.closer_vertices.clear();
        if (on_boundary) {
                for (i=1; i < v->degree-1; i++) {
                        v->vc.farther_vertices.push_back(i);
                }  
                return;
        }
        for (i=0; i < v->degree; i++) {
                calc_next_vertex_geodesic(NULL, v, i, g_next, -1, &g); // only interested in geodesic length
                if (g_next.length > g.length * fudge_factor) {
                        v->vc.farther_vertices.push_back(i);
                } else if (g_next.length < g.length / fudge_factor) {
                        v->vc.closer_vertices.push_back((vertex*)(v->incident_vertices[i]));
                }
                // TODO: Remove
                else {
                        int j=1;
                }
        }

}

/* Given a shortest path from boundary b to vertex v, function
 * calculates the geodesic from boundary b to the vertex
 * v->incident_vertices[k].
 */
void calc_next_vertex_geodesic(mesh* m, vertex* v, int k, geodesic& g, int b, geodesic const* g_source)
{
        int i;
        vertex *v_next;
        edge* e = (edge*)v->incident_edges[k];
        double alpha, beta, offset;
        geodesic path;
        double angle, length, position;
        double sign = 1;
        g.boundary = b;

        if (g_source) {
                path = *g_source;
        }
        else {
                path  = v->shortest_paths[b];
        }
        alpha = calc_geodesic_edge_angle(v, path, k);
        if (alpha > M_PI) {
                sign = -1;
                alpha = 2 * M_PI - alpha;
        }
        g.length = asinh(sinh(path.length) * e->cosh_length - cosh(path.length) *
                       e->sinh_length * cos(alpha) );
        offset = acosh((sinh(g.length) * sinh(path.length) + e->cosh_length ) /
                        (cosh(g.length) * cosh(path.length)));
        beta = acos(-1 * cos(alpha) * cosh(offset) + sin(alpha) * sinh(offset) * sinh(path.length));
        offset *= sign;
        beta *= sign;
        if (m) {
                g.position = normalize_position(m, path.position + offset, b);
        }
        g.angle = calc_next_geodesic_edge_angle(v, (vertex*)(v->incident_vertices[k]), beta);
        g.originating_vertex = v;
}

/* Assuming that geodesic path hits vertex v, function calculuates the
 * angle (0 <= angle < 2PI) of rotation from the geodesic to the k-th
 * outgoing edge incident to v.
 */
double calc_geodesic_edge_angle(vertex *v, const geodesic& path, int k)
{
        double alpha=0;
        for (int i=0; i<k; i++) {
                alpha += *(v->inner_angles[i]);
        }
        alpha += path.angle;
        return normalize_angle(alpha);
}

/* Assuming that the geodesic hitting vertex v_next makes an angle beta
 * with the edge joining v to v_next (measured counterclockwise from
 * the edge to the geodesic), calculates the standard angle for the geodesic -
 * i.e. counterclockwise rotation from the geodesic to the 0th-incident edge.
 */
double calc_next_geodesic_edge_angle(vertex *v, vertex* v_next, double beta)
{
        int i; 
        int k=-1;
        double angle = -1 * beta;
        for (i=0; i<v_next->degree; i++) {
                if (v == (vertex*)(v_next->incident_vertices[i])) {
                        k = i;
                        break;
                }
        }
        if (k < 0) {
                printf("Error finding incident vertices in calc_next_geodesic_edge_angle.\n");
                exit(1);
        }
        for (i=k; i<v_next->degree - v_next->boundary; i++) {
                angle += *(v_next->inner_angles[i]);
        }
        angle += v_next->boundary * M_PI;
        return normalize_angle(angle);
}

/* Adds geodesic g to the list of geodesics from the boundary to v.
 * Makes note of the origin of the geodesic (i.e. which vertex came before)
 * and removes from the list of vertices where geodesics are expected to
 * come from (all closer vertices).
 *
 * If the vertex hasn't been hit with a geodesic before, it uses the
 * geodesic g to determie which incident vertices are closer and farther.
 */
void add_geodesic_to_vertex(vertex* v, vertex* orig_v, const geodesic& g)
{
        if (v->geodesics.empty()) {
                calc_vertex_config(v, g, false);
        }
        v->geodesics.push_back(g);
        v->vc.closer_vertices.remove(orig_v);
}

/* Once the list of geodesics hitting vertex v is generated (from all
 * closer incident vertices) this function calculates their average.
 * In theory they should all be equal, so this function also measures
 * the error between the different geodesics, and returns the result.
 */
geodesic calc_average_geodesic(vertex* v)
{
        int n=0;
        geodesic g, error;
        g = average(v->geodesics);
        error = std_dev(v->geodesics);
        v->shortest_paths[g.boundary] = g;
        return error;
}

/* Given angle, function returns the equivalent angle in 
 * the range 0 <= angle < 2*Pi
 */
double normalize_angle(double angle)
{
        if (angle < 0) {
                return normalize_angle(angle + 2 * M_PI);
        }
        if (angle >= 2 * M_PI) {
                return normalize_angle(angle - 2 * M_PI);
        }
        return angle;
}

/* Given position, function returns the number equivalent
 * to position modulo the length of boundary b. i.e.
 *  0 <= position < boundary length.
 */
double normalize_position(mesh* m, double position, int boundary)
{
        if (position >= m->boundary_lengths[boundary]) {
                return normalize_position(m, position - m->boundary_lengths[boundary], boundary);
        } else if (position < 0) {
                return normalize_position(m, position + m->boundary_lengths[boundary], boundary);
        } else {
                return position;
        }
}

