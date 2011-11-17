CC = g++
DEBUG = -ggdb
CFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	-IIML++ -ISparseLib++/include -ISparseLib++/mv/include \
	-I/usr/local/include
LFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	SparseLib++/lib/libsparse.a SparseLib++/lib/libspblas.a SparseLib++/lib/libmv.a \
	-L/usr/local/lib -llbfgs -lgsl -lgslcblas

HEADERS = \
	   cricci.h		\
	   partition.h		\
	   mesh_improver.h	\
	   ccirclepack.h	\
	   cfile_io.h		\
	   csimpson.h		\
	   geodesic.h		\
	   graph.h		\
	   cmesh.h

OBJS =  \
	   cmain.o		\
	   partition.o		\
	   cricci.o		\
	   mesh_improver.o	\
	   ccirclepack.o	\
	   cfile_io.o		\
	   csimpson.o		\
	   graph.o		\
	   geodesic.o		\
	   cmesh.o

main : $(OBJS)
	$(CC) $(DEBUG) $(LFLAGS) $(OBJS) -o main

cmain.o : cmain.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmain.c

partition.o : partition.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c partition.c

cricci.o : cricci.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cricci.c

mesh_improver.o : mesh_improver.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c mesh_improver.c

ccirclepack.o : ccirclepack.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c ccirclepack.c

cfile_io.o : cfile_io.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cfile_io.c

csimpson.o : csimpson.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c csimpson.c

graph.o : graph.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c graph.c

geodesic.o : geodesic.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c geodesic.c

cmesh.o : cmesh.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmesh.c
