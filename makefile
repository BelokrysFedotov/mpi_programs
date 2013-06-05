DIRS = Test Sum HeatEquation Integral_sin(1_x)

.PHONY : clean build-all

build-all: 
	for dir in $(DIRS); do cd $$dir;make;cd ..; done
    
clean:
	rm -f ~* *.o *.tmp
	for dir in $(DIRS); do cd $$dir;make clean;cd ..; done