CC=g++ -std=c++11
#CODE_SOURCES = .*cpp
#SOURCES= $(CODE_SOURCES)
#LAPACK= -L /cm/shared/apps/lapack/open64/64/3.5.0/ -L /usr/lib -L /cm/shared/apps/blas/gcc/current/lib64/ -llapack -lblas -lm -lgomp -fopenmp -ffast-math -lpthread
#BOOST= -I /cm/shared/apps/boost/1.55.0/include/boost/ -I /cm/shared/apps/boost/1.55.0/include/

#EXECUTABLE=mcfile

#CODE_OBJECTS=\
	main.o	\
	search_for.o	\

#OBJECTS= $(CODE_OBJECTS)

#headers1=search_for.h 

mcfile: main.o search_for.o
	$(CC) -o mcfile main.o search_for.o

main.o : main.cpp search_for.h
	$(CC)  main.cpp -c
search_for.o :search_for.cpp search_for.h
	$(CC) search_for.cpp -c
clean:
	rm -f $(OBJECTS) $(PROG)
	rm -f $(EXECUTABLE)

depend:
	makedepend $(INCLUDE) $(GSL_INC)  $(BLAS_INC) $(GSL_LIB) $(BLAS_LIB) -- -o $(CFLAGS) -- $(SOURCES) 



# DO NOT DELETE

