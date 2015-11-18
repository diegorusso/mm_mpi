include ../../Makefile.mm.inc

VPATH = ../shared

OBJECTS = mm.o check.o prthead.o gendat.o prtspeed.o mxm.o \
          state.o cclock.o reset.o

EXEC    = x.mm

${EXEC}: ${OBJECTS}
	@echo
	${LD} ${LDFLAGS} -o ${EXEC} ${OBJECTS} ${LIBS}

.SUFFIXES: .c .o
.c.o    :
	@echo
	${CC} ${CFLAGS} ${INCS} -c -o $@ $<

clean   :
	/bin/rm -f ${OBJECTS}

veryclean:
	/bin/rm -f ${OBJECTS} ${EXEC}
