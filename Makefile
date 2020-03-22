CXX=		g++
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall #-fsanitize=address -fno-omit-frame-pointer#-Winline 
CPPFLAGS=
INCLUDES=
OBJS=		Output.o CommandLines.o Process_Read.o Assembly.o kmer.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o #ksw2_extz2_sse.o
EXE=		hifiasm
LIBS=		-lz -lpthread -lm #-fsanitize=address -fno-omit-frame-pointer

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

Assembly.o: Assembly.h Process_Read.h kseq.h Overlaps.h kvec.h kdq.h
Assembly.o: CommandLines.h kmer.h Hash_Table.h khash.h POA.h Correct.h
Assembly.o: Levenshtein_distance.h Output.h Trio.h
CommandLines.o: CommandLines.h ketopt.h
Correct.o: Correct.h Hash_Table.h khash.h kmer.h Process_Read.h kseq.h
Correct.o: Overlaps.h kvec.h kdq.h CommandLines.h Levenshtein_distance.h
Correct.o: POA.h Assembly.h #ksw2.h
Hash_Table.o: Hash_Table.h khash.h kmer.h Process_Read.h kseq.h Overlaps.h
Hash_Table.o: kvec.h kdq.h CommandLines.h Correct.h Levenshtein_distance.h
Hash_Table.o: POA.h ksort.h
Levenshtein_distance.o: Levenshtein_distance.h
Output.o: Output.h CommandLines.h
Overlaps.o: Overlaps.h kvec.h kdq.h ksort.h Process_Read.h kseq.h
Overlaps.o: CommandLines.h
POA.o: POA.h Hash_Table.h khash.h kmer.h Process_Read.h kseq.h Overlaps.h
POA.o: kvec.h kdq.h CommandLines.h Correct.h Levenshtein_distance.h
Process_Read.o: Process_Read.h kseq.h Overlaps.h kvec.h kdq.h CommandLines.h
kmer.o: kmer.h Process_Read.h kseq.h Overlaps.h kvec.h kdq.h CommandLines.h
main.o: CommandLines.h Process_Read.h kseq.h Overlaps.h kvec.h kdq.h
main.o: Assembly.h Levenshtein_distance.h
Trio.o: Trio.h khashl.h kthread.h Process_Read.h CommandLines.h
kthread.o: kthread.h
#ksw2_extz2_sse.o: ksw2.h