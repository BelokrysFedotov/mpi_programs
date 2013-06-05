integral: integral.c
	mpicc integral.c -o integral
run: integral
	mpirun -n 1 integral 1.0 10.0 0.01
clean:
	rm -f *.o *.tmp ~* integral