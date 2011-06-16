CC = g++
DEBUG = -ggdb
CFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	-IIML++ -ISparseLib++/include -ISparseLib++/mv/include
LFLAGS = -ansi -pedantic  '-DCOMPLEX=std::complex<double>' \
	SparseLib++/lib/libsparse.a SparseLib++/lib/libspblas.a SparseLib++/lib/libmv.a \
	-L/usr/local/lib -llbfgs

HEADERS = \
	   cricci.h		\
	   ccirclepack.h	\
	   cfile_io.h		\
		 geodesic.h
	   csimpson.h		\
	   cmesh.h

OBJS =  \
	   cmain.o		\
	   cricci.o		\
	   ccirclepack.o	\
	   cfile_io.o		\
	   csimpson.o		\
	   cmesh.o

main : $(OBJS)
	$(CC) $(DEBUG) $(LFLAGS) $(OBJS) -o main -lm

cmain.o : cmain.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmain.c

cricci.o : cricci.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cricci.c

ccirclepack.o : ccirclepack.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c ccirclepack.c

cfile_io.o : cfile_io.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cfile_io.c

csimpson.o : csimpson.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c csimpson.c

cmesh.o : cmesh.c $(HEADERS)
	$(CC) $(DEBUG) $(CFLAGS) -c cmesh.c
