include ../../makeconf

INCPATH = -I../../include -I../THIRD_PARTY
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  binTree.cpp\
		dijkstraGraph.cpp\
		heap.cpp\
		list.cpp\
		clusterGraph.cpp\
		graph.cpp\
		jqsort.cpp\
		matrix.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/primitives.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/primitives.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/primitives.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

binTree.o: /usr/include/stdio.h /usr/include/features.h
binTree.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
binTree.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
binTree.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
binTree.o: /usr/include/bits/typesizes.h /usr/include/libio.h
binTree.o: /usr/include/_G_config.h /usr/include/wchar.h
binTree.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
binTree.o: ../../include/binTree.h
dijkstraGraph.o: ../../include/dijkstraGraph.h ../../include/heap.h
dijkstraGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
dijkstraGraph.o: /usr/include/features.h /usr/include/bits/predefs.h
dijkstraGraph.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
dijkstraGraph.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
dijkstraGraph.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
dijkstraGraph.o: /usr/include/libio.h /usr/include/_G_config.h
dijkstraGraph.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
dijkstraGraph.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
dijkstraGraph.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
dijkstraGraph.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
dijkstraGraph.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
dijkstraGraph.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
dijkstraGraph.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
dijkstraGraph.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
dijkstraGraph.o: ../../include/graph.h ../../include/list.h
heap.o: /usr/include/stdio.h /usr/include/features.h
heap.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
heap.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
heap.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
heap.o: /usr/include/bits/typesizes.h /usr/include/libio.h
heap.o: /usr/include/_G_config.h /usr/include/wchar.h
heap.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
heap.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
heap.o: /usr/include/endian.h /usr/include/bits/endian.h
heap.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
heap.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
heap.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
heap.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
heap.o: ../../include/heap.h ../../include/j_mesh.h /usr/include/math.h
heap.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
heap.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
heap.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
heap.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
heap.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
heap.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
list.o: /usr/include/stdio.h /usr/include/features.h
list.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
list.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
list.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
list.o: /usr/include/bits/typesizes.h /usr/include/libio.h
list.o: /usr/include/_G_config.h /usr/include/wchar.h
list.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
list.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
list.o: /usr/include/endian.h /usr/include/bits/endian.h
list.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
list.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
list.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
list.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
list.o: ../../include/list.h ../../include/jqsort.h
clusterGraph.o: ../../include/clusterGraph.h ../../include/heap.h
clusterGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
clusterGraph.o: /usr/include/features.h /usr/include/bits/predefs.h
clusterGraph.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
clusterGraph.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
clusterGraph.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
clusterGraph.o: /usr/include/libio.h /usr/include/_G_config.h
clusterGraph.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
clusterGraph.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
clusterGraph.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
clusterGraph.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
clusterGraph.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
clusterGraph.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
clusterGraph.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
clusterGraph.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
clusterGraph.o: ../../include/graph.h ../../include/list.h
graph.o: ../../include/graph.h ../../include/list.h /usr/include/stdio.h
graph.o: /usr/include/features.h /usr/include/bits/predefs.h
graph.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
graph.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
graph.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
graph.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
graph.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
matrix.o: /usr/include/math.h /usr/include/features.h
matrix.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
matrix.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
matrix.o: /usr/include/gnu/stubs-64.h /usr/include/bits/huge_val.h
matrix.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
matrix.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
matrix.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
matrix.o: ../../include/matrix.h /usr/include/stdio.h
matrix.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
matrix.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
matrix.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
matrix.o: ../../include/list.h
