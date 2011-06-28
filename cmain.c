// cmain.c
//
#include <stdlib.h>
#include <stdio.h>
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
#include "cricci.h"
#include "graph.h"

int main(int argc, char *argv[]) 
{
        filedata data;
        mesh m;
        if ( argc < 2 ) {
                printf( "usage: %s in_filename [out_filename]\n", argv[0] );
                return 0;
        }
        initialize_filedata(&data);
        read(argv[1], &data);
        initialize_mesh(&m, &data);
        deallocate_filedata(&data);
        calc_circlepack_metric(&m);
        run_ricci_flow(&m);
        if (argc == 3) {
                save_mesh(argv[2], (void *) &m);
        }
        calc_boundary_lengths(&m);
        for (int i=0; i<m.boundary_count; ++i) {
                printf("B%i: %f; ", i, m.boundary_lengths[i]);
        }
        printf("\n");
        calc_vertex_boundary_distances(&m);
        deallocate_mesh(&m);
}
