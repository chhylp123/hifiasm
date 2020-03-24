CXX=		g++
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CPPFLAGS=
INCLUDES=
OBJS=		Output.o CommandLines.o Process_Read.o Assembly.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o
EXE=		hifiasm
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .o
.PHONY:all clean depend

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

Assembly.o: Assembly.h CommandLines.h Process_Read.h kseq.h Overlaps.h kvec.h
Assembly.o: kdq.h kmer.h Hash_Table.h khashl.h POA.h Correct.h
Assembly.o: Levenshtein_distance.h Output.h Trio.h
CommandLines.o: CommandLines.h ketopt.h
Correct.o: Correct.h Hash_Table.h khashl.h kmer.h Process_Read.h kseq.h
Correct.o: Overlaps.h kvec.h kdq.h CommandLines.h Levenshtein_distance.h
Correct.o: POA.h Assembly.h
Hash_Table.o: Hash_Table.h khashl.h kmer.h Process_Read.h kseq.h Overlaps.h
Hash_Table.o: kvec.h kdq.h CommandLines.h Correct.h Levenshtein_distance.h
Hash_Table.o: POA.h ksort.h
Levenshtein_distance.o: Levenshtein_distance.h
Output.o: Output.h CommandLines.h
Overlaps.o: Overlaps.h kvec.h kdq.h ksort.h Process_Read.h kseq.h
Overlaps.o: CommandLines.h Hash_Table.h khashl.h kmer.h Correct.h
Overlaps.o: Levenshtein_distance.h POA.h
POA.o: POA.h Hash_Table.h khashl.h kmer.h Process_Read.h kseq.h Overlaps.h
POA.o: kvec.h kdq.h CommandLines.h Correct.h Levenshtein_distance.h
Process_Read.o: Process_Read.h kseq.h Overlaps.h kvec.h kdq.h CommandLines.h
Trio.o: khashl.h kthread.h Process_Read.h kseq.h Overlaps.h kvec.h kdq.h
Trio.o: CommandLines.h Trio.h kmer.h
kthread.o: kthread.h
main.o: CommandLines.h Process_Read.h kseq.h Overlaps.h kvec.h kdq.h
main.o: Assembly.h Levenshtein_distance.h
