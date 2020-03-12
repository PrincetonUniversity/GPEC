build:
	$(F90) -ldl dll_module.f90 cvode.f90 main.f90 -c
	$(F90) -o a.out cvode.o dll_module.o main.o -ldl
	gcc -c -Wall -fpic xintgrnd.c
	gcc -shared -o libxintgrnd.so xintgrnd.o
	mkdir -p ../lib
	ar -r ../lib/libcvode.a *.o
	mkdir -p ../bin

clean:
	rm -f *.o *.mod *.out
