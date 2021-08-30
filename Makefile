CXX=		g++
CC=			gcc
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		CommandLines.o Process_Read.o Assembly.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o Purge_Dups.o \
			htab.o hist.o sketch.o anchor.o extract.o sys.o ksw2_extz2_sse.o hic.o rcut.o horder.o \
			tovlp.o inter.o
EXE=		hifiasm
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

Assembly.o: Assembly.h CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h
Assembly.o: Hash_Table.h htab.h POA.h Correct.h Levenshtein_distance.h
Assembly.o: kthread.h
CommandLines.o: CommandLines.h ketopt.h
Correct.o: Correct.h Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h
Correct.o: kdq.h CommandLines.h Levenshtein_distance.h POA.h Assembly.h
Correct.o: ksw2.h ksort.h
Hash_Table.o: Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h kdq.h
Hash_Table.o: CommandLines.h ksort.h
Levenshtein_distance.o: Levenshtein_distance.h
Output.o: Output.h CommandLines.h
Overlaps.o: Overlaps.h kvec.h kdq.h ksort.h Process_Read.h CommandLines.h
Overlaps.o: Hash_Table.h htab.h Correct.h Levenshtein_distance.h POA.h
Overlaps.o: Purge_Dups.h
POA.o: POA.h Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h kdq.h
POA.o: CommandLines.h Correct.h Levenshtein_distance.h
Process_Read.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
Purge_Dups.o: ksort.h Purge_Dups.h kvec.h kdq.h Overlaps.h Hash_Table.h
Purge_Dups.o: htab.h Process_Read.h CommandLines.h Correct.h
Purge_Dups.o: Levenshtein_distance.h POA.h kthread.h
Trio.o: khashl.h kthread.h kseq.h Process_Read.h Overlaps.h kvec.h kdq.h
Trio.o: CommandLines.h htab.h
anchor.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
anchor.o: ksort.h Hash_Table.h
extract.o: Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h khashl.h
extract.o: kseq.h
hist.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
htab.o: kthread.h khashl.h kseq.h ksort.h htab.h Process_Read.h Overlaps.h
htab.o: kvec.h kdq.h CommandLines.h
kthread.o: kthread.h
main.o: CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h Assembly.h
main.o: Levenshtein_distance.h htab.h
sketch.o: kvec.h htab.h Process_Read.h Overlaps.h kdq.h CommandLines.h
sys.o: htab.h Process_Read.h Overlaps.h kvec.h kdq.h CommandLines.h
hic.o: hic.h
rcut.o: rcut.h
horder.o: horder.h
tovlp.o: tovlp.h
inter.o: inter.h
