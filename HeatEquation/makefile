heateq: heateq.c
	mpicc heateq.c -o heateq
run: heateq
	./mpirun -n 1 heateq
clean:
	rm -f *.o *.tmp ~* heateq