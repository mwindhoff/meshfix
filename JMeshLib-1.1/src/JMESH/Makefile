include ../../makeconf

INCPATH = -I../../include

CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  jmesh.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/jmesh.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/jmesh.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/jmesh.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

jmesh.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
jmesh.o: /usr/include/features.h /usr/include/bits/predefs.h
jmesh.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
jmesh.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
jmesh.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
jmesh.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
jmesh.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
jmesh.o: /usr/include/math.h /usr/include/bits/huge_val.h
jmesh.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
jmesh.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
jmesh.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
jmesh.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
jmesh.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
jmesh.o: /usr/include/bits/posix2_lim.h ../../include/list.h
jmesh.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
jmesh.o: ../../include/triangle.h ../../include/matrix.h ../../include/tin.h
jmesh.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
jmesh.o: /usr/include/endian.h /usr/include/bits/endian.h
jmesh.o: /usr/include/bits/byteswap.h /usr/include/sys/select.h
jmesh.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
jmesh.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
jmesh.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
jmesh.o: /usr/include/string.h /usr/include/xlocale.h
