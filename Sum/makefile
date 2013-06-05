sum: sum.c
	mpicc sum.c -o sum
run: sum
	./mpirun -n 1 sum
clean:
	rm -f *.o *.tmp ~* sum