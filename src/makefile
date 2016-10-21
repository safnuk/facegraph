CC = g++
DEBUG = -ggdb
SPARSE_LIB_DIR = SparseLib++

CFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	-IIML++ -I$(SPARSE_LIB_DIR)/include -I$(SPARSE_LIB_DIR)/mv/include \
	-I/usr/local/include
LFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	$(SPARSE_LIB_DIR)/lib/libsparse.a $(SPARSE_LIB_DIR)/lib/libspblas.a $(SPARSE_LIB_DIR)/lib/libmv.a \
	-llbfgs -lgsl -lgslcblas -lm


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

main : $(OBJS) sparselib
	$(CC) $(DEBUG) $(LFLAGS) $(OBJS) -o main

sparselib :
	$(MAKE) -C $(SPARSE_LIB_DIR) sp

cmain.o : cmain.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmain.cc

partition.o : partition.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c partition.cc

cricci.o : cricci.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cricci.cc

mesh_improver.o : mesh_improver.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c mesh_improver.cc

ccirclepack.o : ccirclepack.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c ccirclepack.cc

cfile_io.o : cfile_io.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cfile_io.cc

csimpson.o : csimpson.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c csimpson.cc

graph.o : graph.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c graph.cc

geodesic.o : geodesic.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c geodesic.cc

cmesh.o : cmesh.cc $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmesh.cc
