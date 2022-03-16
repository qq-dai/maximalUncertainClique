# all: uclique
# .PHONY : all

CC = g++ "-std=c++11"
CONFIG = -O3 -g

objects = algorithms.o maximalClique.o main.o

.PHONY : clean

uclique: $(objects)
	$(CC) -o uclique $(objects) $(CONFIG)

%.o:%.cpp
	$(CC) -c $^ $(CONFIG)

clean:
	rm -f *.o uclique