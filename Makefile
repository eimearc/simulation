CC = g++
CFLAGS = -std=c++17
LDFLAGS = `pkg-config --static --libs glfw3`

sim: main.cpp
	$(CC) -o $@ main.cpp $(CFLAGS) $(LDFLAGS)

clean:
	rm sim
