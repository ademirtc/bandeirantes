all: main runexp

#Compiladores
CC=gcc
CXX=g++

FLAGS= -O3 -Wall

#Bibliotecas
GFTLIB  = -L$(GFT_DIR)/lib -lgft
GFTFLAGS  = -I$(GFT_DIR)/include

MUNKLIB   = -L./lib/munkres-cpp/build -lmunkres
MUNKFLAGS = -I./lib/munkres-cpp/src -std=gnu++11


#Rules
libgft:
	$(MAKE) -C $(GFT_DIR)

main: main.cpp libgft
	$(CXX) $(FLAGS) $(GFTFLAGS) $(MUNKFLAGS) \
	main.cpp $(GFTLIB) $(MUNKLIB) -o main -lm -lz 

runexp: runexp.cpp libgft
	$(CXX) $(FLAGS) $(GFTFLAGS) \
	runexp.cpp $(GFTLIB) -o runexp -lm -lz 

clean:
	$(RM) *~ *.o main runexp


