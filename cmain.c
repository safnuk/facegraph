// cmain.c
//
#include <stdlib.h>
#include <stdio.h>
#include "cfile_io.h"
#include "cmesh.h"
#include "ccirclepack.h"

int main(int argc, char *argv[]) {
        filedata data;
        mesh m;
        if ( argc != 2 ) {
                printf( "usage: %s filename\n", argv[0] );
                return 0;
        }
        initialize_filedata(&data);
        read(argv[1], &data);
        initialize_mesh(&m, &data);
        deallocate_filedata(&data);
        print_mesh(&m);
        calc_circlepack_metric(&m);
        deallocate_mesh(&m);
}
