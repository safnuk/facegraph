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
#include "cmesh.h"
#include "ccirclepack.h"
#include "cricci.h"

int main(int argc, char *argv[]) 
{
        filedata data;
        mesh m;
        mesh m_double;
        if ( argc < 2 ) {
                printf( "usage: %s in_filename [out_filename]\n", argv[0] );
                return 0;
        }
        initialize_filedata(&data);
        read(argv[1], &data);
        initialize_mesh(&m, &data);
        deallocate_filedata(&data);
        double_mesh(&m, &m_double);
        calc_circlepack_metric(&m_double);
        run_ricci_flow(&m_double);
        split_doubled_mesh(&m, &m_double);
        if (argc == 3) {
                save_mesh(argv[2], (void *) &m);
        }
        deallocate_mesh(&m);
        deallocate_mesh(&m_double);
}
