CC = g++
CFLAGS = -std=c++17
GLEW_PATH = /usr/local/lib/
LDFLAGS = `pkg-config --static --libs glfw3` -L$(GLEW_PATH) -lGLEW -framework OpenGL

sim: main.cpp
	$(CC) -o $@ main.cpp $(CFLAGS) $(LDFLAGS)

clean:
	rm sim
