BUILD_DIR?= build
SRC_DIR?=   src
CXX=		g++
CC=			gcc
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		$(BUILD_DIR)/CommandLines.o $(BUILD_DIR)/Process_Read.o $(BUILD_DIR)/Assembly.o $(BUILD_DIR)/Hash_Table.o \
			$(BUILD_DIR)/POA.o $(BUILD_DIR)/Correct.o $(BUILD_DIR)/Levenshtein_distance.o $(BUILD_DIR)/Overlaps.o $(BUILD_DIR)/Trio.o $(BUILD_DIR)/kthread.o $(BUILD_DIR)/Purge_Dups.o \
			$(BUILD_DIR)/htab.o $(BUILD_DIR)/hist.o $(BUILD_DIR)/sketch.o $(BUILD_DIR)/anchor.o $(BUILD_DIR)/extract.o $(BUILD_DIR)/sys.o $(BUILD_DIR)/hic.o $(BUILD_DIR)/rcut.o $(BUILD_DIR)/horder.o \
			$(BUILD_DIR)/tovlp.o $(BUILD_DIR)/inter.o $(BUILD_DIR)/kalloc.o $(BUILD_DIR)/gfa_ut.o $(BUILD_DIR)/gchain_map.o
EXE=		$(BUILD_DIR)/hifiasm
LIBS=		-lz -lpthread -lm

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
.PHONY:all clean depend

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
		mkdir -p $(BUILD_DIR); $(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
		mkdir -p $(BUILD_DIR); $(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(EXE)

$(EXE):$(OBJS) $(BUILD_DIR)/main.o
		mkdir -p $(BUILD_DIR); $(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -rf $(BUILD_DIR) gmon.out $(EXE) *~ *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

$(BUILD_DIR)/Assembly.o: $(SRC_DIR)/Assembly.h $(SRC_DIR)/CommandLines.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h
$(BUILD_DIR)/Assembly.o: $(SRC_DIR)/Hash_Table.h $(SRC_DIR)/htab.h $(SRC_DIR)/POA.h $(SRC_DIR)/Correct.h $(SRC_DIR)/Levenshtein_distance.h
$(BUILD_DIR)/Assembly.o: $(SRC_DIR)/kthread.h
$(BUILD_DIR)/CommandLines.o: $(SRC_DIR)/CommandLines.h $(SRC_DIR)/ketopt.h
$(BUILD_DIR)/Correct.o: $(SRC_DIR)/Correct.h $(SRC_DIR)/Hash_Table.h $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h
$(BUILD_DIR)/Correct.o: $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h $(SRC_DIR)/Levenshtein_distance.h $(SRC_DIR)/POA.h $(SRC_DIR)/Assembly.h
$(BUILD_DIR)/Correct.o: $(SRC_DIR)/ksort.h
$(BUILD_DIR)/Hash_Table.o: $(SRC_DIR)/Hash_Table.h $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h
$(BUILD_DIR)/Hash_Table.o: $(SRC_DIR)/CommandLines.h $(SRC_DIR)/ksort.h
$(BUILD_DIR)/Levenshtein_distance.o: $(SRC_DIR)/Levenshtein_distance.h
$(BUILD_DIR)/Output.o: $(SRC_DIR)/Output.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/Overlaps.o: $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/ksort.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/Overlaps.o: $(SRC_DIR)/Hash_Table.h $(SRC_DIR)/htab.h $(SRC_DIR)/Correct.h $(SRC_DIR)/Levenshtein_distance.h $(SRC_DIR)/POA.h
$(BUILD_DIR)/Overlaps.o: $(SRC_DIR)/Purge_Dups.h
$(BUILD_DIR)/POA.o: $(SRC_DIR)/POA.h $(SRC_DIR)/Hash_Table.h $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h
$(BUILD_DIR)/POA.o: $(SRC_DIR)/CommandLines.h $(SRC_DIR)/Correct.h $(SRC_DIR)/Levenshtein_distance.h
$(BUILD_DIR)/Process_Read.o: $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/Purge_Dups.o: $(SRC_DIR)/ksort.h $(SRC_DIR)/Purge_Dups.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/Hash_Table.h
$(BUILD_DIR)/Purge_Dups.o: $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/CommandLines.h $(SRC_DIR)/Correct.h
$(BUILD_DIR)/Purge_Dups.o: $(SRC_DIR)/Levenshtein_distance.h $(SRC_DIR)/POA.h $(SRC_DIR)/kthread.h
$(BUILD_DIR)/Trio.o: $(SRC_DIR)/khashl.h $(SRC_DIR)/kthread.h $(SRC_DIR)/kseq.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h
$(BUILD_DIR)/Trio.o: $(SRC_DIR)/CommandLines.h $(SRC_DIR)/htab.h
$(BUILD_DIR)/anchor.o: $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/anchor.o: $(SRC_DIR)/ksort.h $(SRC_DIR)/Hash_Table.h
$(BUILD_DIR)/extract.o: $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h $(SRC_DIR)/khashl.h
$(BUILD_DIR)/extract.o: $(SRC_DIR)/kseq.h
$(BUILD_DIR)/hist.o: $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/htab.o: $(SRC_DIR)/kthread.h $(SRC_DIR)/khashl.h $(SRC_DIR)/kseq.h $(SRC_DIR)/ksort.h $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h
$(BUILD_DIR)/htab.o: $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/kthread.o: $(SRC_DIR)/kthread.h
$(BUILD_DIR)/main.o: $(SRC_DIR)/CommandLines.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/Assembly.h
$(BUILD_DIR)/main.o: $(SRC_DIR)/Levenshtein_distance.h $(SRC_DIR)/htab.h
$(BUILD_DIR)/sketch.o: $(SRC_DIR)/kvec.h $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/sys.o: $(SRC_DIR)/htab.h $(SRC_DIR)/Process_Read.h $(SRC_DIR)/Overlaps.h $(SRC_DIR)/kvec.h $(SRC_DIR)/kdq.h $(SRC_DIR)/CommandLines.h
$(BUILD_DIR)/hic.o: $(SRC_DIR)/hic.h
$(BUILD_DIR)/rcut.o: $(SRC_DIR)/rcut.h
$(BUILD_DIR)/horder.o: $(SRC_DIR)/horder.h
$(BUILD_DIR)/tovlp.o: $(SRC_DIR)/tovlp.h
$(BUILD_DIR)/inter.o: $(SRC_DIR)/inter.h $(SRC_DIR)/Process_Read.h
$(BUILD_DIR)/kalloc.o: $(SRC_DIR)/kalloc.h
$(BUILD_DIR)/gfa_ut.o: $(SRC_DIR)/Overlaps.h
$(BUILD_DIR)/gchain_map.o: $(SRC_DIR)/gchain_map.h
