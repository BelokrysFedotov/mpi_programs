test: test.c
	mpicc test.c -o test
run: test
	./mpirun -n 1 test
clean:
	rm -f *.o *.tmp ~* test