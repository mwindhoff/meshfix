include ../../makeconf

INCPATH = -I../../include
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH) -I../RANGECODER

SOURCES.C    =  point.cpp\
		vertex.cpp\
		edge.cpp\
		triangle.cpp\
		tin.cpp\
		io.cpp\
		checkAndRepair.cpp
		
OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/mesh_structure.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/mesh_structure.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/mesh_structure.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

point.o: ../../include/point.h ../../include/j_mesh.h /usr/include/stdio.h
point.o: /usr/include/features.h /usr/include/bits/predefs.h
point.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
point.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
point.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
point.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
point.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
point.o: /usr/include/math.h /usr/include/bits/huge_val.h
point.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
point.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
point.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
point.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
point.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
point.o: /usr/include/bits/posix2_lim.h /usr/include/stdlib.h
point.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/endian.h
point.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
point.o: /usr/include/sys/select.h /usr/include/bits/select.h
point.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
point.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
point.o: /usr/include/alloca.h /usr/include/errno.h /usr/include/bits/errno.h
point.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
point.o: /usr/include/asm-generic/errno.h
point.o: /usr/include/asm-generic/errno-base.h
vertex.o: ../../include/vertex.h ../../include/j_mesh.h /usr/include/stdio.h
vertex.o: /usr/include/features.h /usr/include/bits/predefs.h
vertex.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
vertex.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
vertex.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
vertex.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
vertex.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
vertex.o: /usr/include/math.h /usr/include/bits/huge_val.h
vertex.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
vertex.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
vertex.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
vertex.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
vertex.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
vertex.o: /usr/include/bits/posix2_lim.h ../../include/list.h
vertex.o: ../../include/point.h ../../include/edge.h ../../include/triangle.h
vertex.o: ../../include/jmesh.h ../../include/matrix.h ../../include/tin.h
vertex.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
vertex.o: /usr/include/endian.h /usr/include/bits/endian.h
vertex.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
vertex.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
vertex.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
vertex.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
vertex.o: /usr/include/errno.h /usr/include/bits/errno.h
vertex.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
vertex.o: /usr/include/asm-generic/errno.h
vertex.o: /usr/include/asm-generic/errno-base.h
edge.o: ../../include/edge.h ../../include/j_mesh.h /usr/include/stdio.h
edge.o: /usr/include/features.h /usr/include/bits/predefs.h
edge.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
edge.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
edge.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
edge.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
edge.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
edge.o: /usr/include/math.h /usr/include/bits/huge_val.h
edge.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
edge.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
edge.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
edge.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
edge.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
edge.o: /usr/include/bits/posix2_lim.h ../../include/vertex.h
edge.o: ../../include/list.h ../../include/point.h ../../include/triangle.h
edge.o: ../../include/jmesh.h ../../include/matrix.h ../../include/tin.h
triangle.o: ../../include/triangle.h ../../include/jmesh.h
triangle.o: ../../include/j_mesh.h /usr/include/stdio.h
triangle.o: /usr/include/features.h /usr/include/bits/predefs.h
triangle.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
triangle.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
triangle.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
triangle.o: /usr/include/libio.h /usr/include/_G_config.h
triangle.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
triangle.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
triangle.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
triangle.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
triangle.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
triangle.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
triangle.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
triangle.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
triangle.o: ../../include/list.h ../../include/point.h ../../include/vertex.h
triangle.o: ../../include/edge.h ../../include/matrix.h ../../include/tin.h
triangle.o: /usr/include/stdlib.h /usr/include/sys/types.h
triangle.o: /usr/include/time.h /usr/include/endian.h
triangle.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
triangle.o: /usr/include/sys/select.h /usr/include/bits/select.h
triangle.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
triangle.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
triangle.o: /usr/include/alloca.h
tin.o: ../../include/tin.h ../../include/jmesh.h ../../include/j_mesh.h
tin.o: /usr/include/stdio.h /usr/include/features.h
tin.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
tin.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
tin.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
tin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
tin.o: /usr/include/_G_config.h /usr/include/wchar.h
tin.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
tin.o: /usr/include/math.h /usr/include/bits/huge_val.h
tin.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
tin.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
tin.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
tin.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
tin.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
tin.o: /usr/include/bits/posix2_lim.h ../../include/list.h
tin.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
tin.o: ../../include/triangle.h ../../include/matrix.h /usr/include/stdlib.h
tin.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/endian.h
tin.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
tin.o: /usr/include/sys/select.h /usr/include/bits/select.h
tin.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
tin.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
tin.o: /usr/include/alloca.h /usr/include/string.h /usr/include/xlocale.h
io.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
io.o: /usr/include/features.h /usr/include/bits/predefs.h
io.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
io.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
io.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
io.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
io.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
io.o: /usr/include/math.h /usr/include/bits/huge_val.h
io.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
io.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
io.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
io.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
io.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
io.o: /usr/include/bits/posix2_lim.h ../../include/list.h
io.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
io.o: ../../include/triangle.h ../../include/matrix.h ../../include/tin.h
io.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
io.o: /usr/include/endian.h /usr/include/bits/endian.h
io.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
io.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
io.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
io.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
io.o: /usr/include/string.h /usr/include/xlocale.h /usr/include/ctype.h
checkAndRepair.o: ../../include/jmesh.h ../../include/j_mesh.h
checkAndRepair.o: /usr/include/stdio.h /usr/include/features.h
checkAndRepair.o: /usr/include/bits/predefs.h /usr/include/sys/cdefs.h
checkAndRepair.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
checkAndRepair.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
checkAndRepair.o: /usr/include/bits/typesizes.h /usr/include/libio.h
checkAndRepair.o: /usr/include/_G_config.h /usr/include/wchar.h
checkAndRepair.o: /usr/include/bits/stdio_lim.h
checkAndRepair.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
checkAndRepair.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
checkAndRepair.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
checkAndRepair.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
checkAndRepair.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
checkAndRepair.o: /usr/include/bits/posix1_lim.h
checkAndRepair.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
checkAndRepair.o: /usr/include/bits/posix2_lim.h ../../include/list.h
checkAndRepair.o: ../../include/point.h ../../include/vertex.h
checkAndRepair.o: ../../include/edge.h ../../include/triangle.h
checkAndRepair.o: ../../include/matrix.h ../../include/tin.h
checkAndRepair.o: ../../include/jqsort.h /usr/include/stdlib.h
checkAndRepair.o: /usr/include/sys/types.h /usr/include/time.h
checkAndRepair.o: /usr/include/endian.h /usr/include/bits/endian.h
checkAndRepair.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
checkAndRepair.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
checkAndRepair.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
checkAndRepair.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
checkAndRepair.o: /usr/include/string.h /usr/include/xlocale.h
