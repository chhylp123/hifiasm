CC=g++

CFLAGS = -w -c -msse4.2 -mpopcnt -fomit-frame-pointer -Winline -O3 -lz
LDFLAGS = -lm -lz -lpthread -O3 -mpopcnt -msse4.2 -lz -w

SOURCES = main.cpp Output.cpp CommandLines.cpp Process_Read.cpp Assembly.cpp kmer.cpp Hash_Table.cpp POA.cpp Correct.cpp Levenshtein_distance.cpp edlib.cpp Overlaps.cpp
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = ccs_assembly


all: $(SOURCES) $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* ccs_assembly