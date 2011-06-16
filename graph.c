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
 */
void calc_vertex_boundary_distances(mesh *m)
{
        std::list<vertex*> active;
        int i;
        for (i=0; i<m->boundary_count; i++) {
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
                set_geodesic(v->shortest_path, b, 0, position, M_PI_2);
                calc_vertex_config(v, v->shortest_path, true;);
                active.push_end(v);
                position += acosh(e->cosh_length);
        }
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
        vertex* v, v1;
        geodesic g;
        while (i!=active.end()) {
                v = *i;
                std::list<int>::iterator j = v->vc.farther_vertices.begin();
                for (;j!=v->vc.farther_vertices.end(); j++) {
                        calc_next_vertex_geodesic(v, *j, g, b);
                        v1 = (vertex*)(v->incident_vertices[*j]);
                        add_geodesic_to_vertex(v1, g);
                        if(v1->vc.closer_vertices.empty()) {
                                calc_average_geodesic(v1);
                                active.push_end(v1);
                        }
                }
                i = active.erase(i);
        }
}

/* Calculates which incident vertices are further away from v 
 * and which are closer. The calculation is simpler for
 * points known to be on the boundary, hence the flag on_boundary. 
 */
void calc_vertex_config(vertex* v, const geodesic &g, bool on_boundary=false)
{
}

void calc_next_vertex_geodesic(vertex *v, int i, geodesic &g, int b)
{
        vertex *v_next;
        edge e;
        double angle, length, position;

}

void add_geodesic_to_vertex(vertex *v, const geodesic &g)
{
}
