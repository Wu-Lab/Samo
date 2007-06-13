CC = g++
LIBS =
CFLAGS = -O3 -Wall

#sources
HEADERS = FibHeap.h  List.h  Matrix.h  MultiAlign.h  PairAlign.h  PDB.h  Samo.h  ProteinChain.h  SVD.h  Tree.h  Utils.h
SRCS = FibHeap.cpp  MultiAlign.cpp  PairAlign.cpp  PDB.cpp  Samo.cpp  ProteinChain.cpp  SVD.cpp  Utils.cpp
LIB = libsamo.a
OBJS = $(SRCS:.cpp=.o)

#application
APPS = samo

all: $(APPS)

$(APPS): %: Main.cpp $(LIB) $(SRCS) $(HEADERS) Makefile
	$(CC) $(CFLAGS) $< $(LIB) $(LIBS) -o $@

$(LIB): $(OBJS)
	$(AR) rus $(LIB) $(OBJS)

$(OBJS): %.o: %.cpp $(HEADERS) Makefile
	$(CC) -c $(CFLAGS) $<

clean:
	rm $(APPS) $(LIB) $(OBJS)
