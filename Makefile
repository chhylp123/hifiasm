CXX=		g++
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CPPFLAGS=
INCLUDES=
OBJS=		Output.o CommandLines.o Process_Read.o Assembly.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o \
			htab.o hist.o sketch.o sys.o
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

Assembly.o: Assembly.h CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h
Assembly.o: kmer.h Hash_Table.h htab.h POA.h Correct.h Levenshtein_distance.h
Assembly.o: Output.h
CommandLines.o: CommandLines.h ketopt.h
Correct.o: Correct.h Hash_Table.h kmer.h Process_Read.h Overlaps.h kvec.h
Correct.o: kdq.h CommandLines.h htab.h Levenshtein_distance.h POA.h
Correct.o: Assembly.h
Hash_Table.o: Hash_Table.h kmer.h Process_Read.h Overlaps.h kvec.h kdq.h
Hash_Table.o: CommandLines.h htab.h Correct.h Levenshtein_distance.h POA.h
Hash_Table.o: ksort.h
Levenshtein_distance.o: Levenshtein_distance.h
Output.o: Output.h CommandLines.h
Overlaps.o: Overlaps.h kvec.h kdq.h ksort.h Process_Read.h CommandLines.h
Overlaps.o: Hash_Table.h kmer.h htab.h Correct.h Levenshtein_distance.h POA.h
POA.o: POA.h Hash_Table.h kmer.h Process_Read.h Overlaps.h kvec.h kdq.h
POA.o: CommandLines.h htab.h Correct.h Levenshtein_distance.h
Process_Read.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
Trio.o: khashl.h kthread.h Process_Read.h Overlaps.h kvec.h kdq.h
Trio.o: CommandLines.h htab.h
hist.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
htab.o: kthread.h khashl.h kseq.h ksort.h htab.h Process_Read.h Overlaps.h
htab.o: kvec.h kdq.h CommandLines.h
kthread.o: kthread.h
main.o: CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h Assembly.h
main.o: Levenshtein_distance.h htab.h
sketch.o: kvec.h htab.h Process_Read.h Overlaps.h kdq.h CommandLines.h
sys.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
