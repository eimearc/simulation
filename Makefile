CC = g++
CFLAGS = -std=c++17

sim: main.cpp
	$(CC) -o $@ main.cpp $(CFLAGS)

clean:
	rm sim
