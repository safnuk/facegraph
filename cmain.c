// cmain.c
//
#include <stdlib.h>
#include <stdio.h>
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
#include "mesh_improver.h"
#include "ccirclepack.h"
#include "cricci.h"
#include "graph.h"

int main(int argc, char *argv[]) 
{
  filedata data;
  mesh m;
  ribbon_graph gamma;
  if ( argc != 2 ) {
    printf( "usage: %s in_filename\n", argv[0] );
    return 0;
  }
  initialize_filedata(&data);
  read(argv[1], &data);
  initialize_mesh(&m, &data);
  calc_circlepack_metric(&m);
  int count = 0;
  while (find_problem_vertices(&m, &data) && (count < 4)) {
    ++count;
  }
  calc_circlepack_metric(&m);  // shouldn't be necessary, but
                               // data is lost in find_problem_vertices
  run_ricci_flow(&m);
  calc_boundary_lengths(&m);
  for (int i=0; i<m.boundary_count; ++i) {
    printf("B%i: %f; ", i, m.boundary_lengths[i]);
  }
  printf("\n");
  calc_cutlocus_graph(&m, &gamma);
  print_ribbon_graph(&gamma);
  save_mesh(argv[1], (void *) &m);

  deallocate_filedata(&data);
  deallocate_mesh(&m);
}
