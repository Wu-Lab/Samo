CC = g++
LIBS =
CFLAGS = -O3 -Wall

#sources
HEADERS = AlignParams.h  FibHeap.h  Matrix.h  MemLeak.h  MultiAlign.h  Options.h  PairAlign.h  PDB.h  ProteinChain.h  Samo.h  SVD.h  Utils.h
SRCS = FibHeap.cpp  MultiAlign.cpp  Options.cpp  PairAlign.cpp  PDB.cpp  ProteinChain.cpp  Samo.cpp  SVD.cpp  Utils.cpp
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
