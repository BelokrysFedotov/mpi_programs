flops: flops.f
	mpif77 flops.f -o flops
run: flops
	./mpirun -n 1 flops
clean:
	rm -f *.o *.tmp ~* flops
