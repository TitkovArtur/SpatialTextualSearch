OS := $(shell uname)
ifeq ($(OS),Darwin)
        #CC      = g++
	CC	= /usr/bin/g++
        CFLAGS  = -c -W -Wall -O2 -mavx -std=c++14 -w -march=native
       	LDFLAGS =
else
        CC      = g++
        CFLAGS  = -O3 -mavx -std=c++14 -w -march=native
        LDFLAGS =
endif
SOURCES=BRQ_main.cpp PageFile.cpp relation.cpp inverted_index.cpp tools.cpp queries.cpp GenericTool.cpp mt19937ar.cpp BRQ_Queries.cpp ResultSet.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: BRQ

BRQ: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o run_BRQ $(LDADD)
	
BkQ: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o run_BkQ.exec $(LDADD)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o
	rm -rf *.exec
