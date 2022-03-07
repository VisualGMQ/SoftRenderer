SRCS 	  = $(wildcard *.cpp)
BINS 	  = $(patsubst %.cpp, %.out, ${SRCS})
SDL_DEP   = `sdl2-config --libs --cflags` `pkg-config --libs --cflags glm`
CXX_STD   = -std=c++17 -O2
CXX_FLAGS = -Wall -Wextra

all:${BINS}


%.out: %.cpp renderer.hpp interactive.hpp
	$(CXX) $< -o $@ ${SDL_DEP} ${CXX_STD} ${CXX_FLAGS}

.PHONY:clean
clean:
	-rm *.out
	-rm *.bmp
	-rm -rf *.dSYM
