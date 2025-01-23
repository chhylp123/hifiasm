CXX=		g++
CC=			gcc
# -g：デバッグ情報の作成を有効にする
# -O3：強力な最適化を行う
# -fomit-frame-pointer：フレーム・ポインタを必要としない関数においては、 フレーム・ポインタをレジスタ内に保持しません。 これにより、 フレーム・ポインタの待避、 セットアップ、 復元を行う命令を使わずに済むようになります。 また、 多くの関数において、 レジスタを余分に利用することができるようになります。 
# -Wall：警告をすべて出力
CXXFLAGS=	-g -O3 -msse4.2 -mpopcnt -fomit-frame-pointer -Wall
CFLAGS=		$(CXXFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		CommandLines.o Process_Read.o Assembly.o Hash_Table.o \
			POA.o Correct.o Levenshtein_distance.o Overlaps.o Trio.o kthread.o Purge_Dups.o \
			htab.o hist.o sketch.o anchor.o extract.o sys.o hic.o rcut.o horder.o ecovlp.o\
			tovlp.o inter.o kalloc.o gfa_ut.o gchain_map.o
EXE=		hifiasm
# ライブラリのリンク。lz：zlib、lpthread：マルチスレッド対応、lm：算術演算
LIBS=		-lz -lpthread -lm

# address sanitizerがオンになっていればコマンドにオプション（-fsanitize=address）を追加する
ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.cpp .c .o
# 疑似ターゲット（.PHONY）を宣言
# もしタスク名と同名のファイルやディレクトリがあると混乱するため
.PHONY:all clean depend

# $@:ターゲットファイル名、$<:最初の依存するファイルの名前
.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# .c.o: というターゲットは，.oというファイルが必要になれば，これを.cからつくる というルールである
.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

# makeコマンドのあとに何も書かない場合、EXE（=hifiasm）を実行
all:$(EXE)

# hifiasmを実行した場合の動作
$(EXE):$(OBJS) main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(EXE) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CPPFLAGS) $(DFLAGS) -- *.cpp)

# DO NOT DELETE

Assembly.o: Assembly.h CommandLines.h Process_Read.h Overlaps.h kvec.h kdq.h
Assembly.o: Hash_Table.h htab.h POA.h Correct.h Levenshtein_distance.h
Assembly.o: kthread.h ecovlp.h
CommandLines.o: CommandLines.h ketopt.h
Correct.o: Correct.h Hash_Table.h htab.h Process_Read.h Overlaps.h kvec.h
Correct.o: kdq.h CommandLines.h Levenshtein_distance.h POA.h Assembly.h
Correct.o: ksort.h
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
ecovlp.o: Hash_Table.h Process_Read.h Overlaps.h kthread.h
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
inter.o: inter.h Process_Read.h
kalloc.o: kalloc.h
gfa_ut.o: Overlaps.h
gchain_map.o: gchain_map.h
