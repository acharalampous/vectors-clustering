CC = g++
CFLAGS = -g3 -Wall
OBJ = main.o utils.o dataset.o clusters.o cl_algorithms.o lsh.o metrics.o
EXE = cluster

$(EXE): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@

main.o: main.cpp
	$(CC) $(CFLAGS) -c $^

utils.o: utils.cpp
	$(CC) $(CFLAGS) -c $^

dataset.o: dataset.cpp
	$(CC) $(CFLAGS) -c $^

clusters.o: clusters.cpp
	$(CC) $(CFLAGS) -c $^

cl_algorithms.o: cl_algorithms.cpp
	$(CC) $(CFLAGS) -c $^

lsh.o: lsh.cpp
	$(CC) $(CFLAGS) -c $^

metrics.o: metrics.cpp
	$(CC) $(CFLAGS) -c $^

.PHONY: clean

clean:
	rm -f $(EXE) $(OBJ) 
