
default: all

serial:
	cd src/serial && make veryclean && make && cp x.mm ../../x.mm_serial;

2D-cannon:
	cd src/cannon && make veryclean && make && cp x.mm ../../x.mm_2D_cannon;

2D-cannon-cblas:
	cd src/cannon && make -f Makefile.cblas veryclean && make -f Makefile.cblas && cp x.mm ../../x.mm_2D_cannon_cblas;

2D-cannon-nonblock:
	cd src/cannon && make -f Makefile.nonblock veryclean && make -f Makefile.nonblock && cp x.mm ../../x.mm_2D_cannon_nonblock;

2D-cannon-nonblock-cblas:
	cd src/cannon && make -f Makefile.nonblock.cblas veryclean && make -f Makefile.nonblock.cblas && cp x.mm ../../x.mm_2D_cannon_nonblock_cblas;

2D-cannon-nonblock-cblas-prepend:
	cd src/cannon && make -f Makefile.nonblock.cblas.prepend veryclean && make -f Makefile.nonblock.cblas.prepend && cp x.mm ../../x.mm_2D_cannon_nonblock_cblas_prepend;

all : serial 2D-cannon 2D-cannon-nonblock 2D-cannon-cblas 2D-cannon-nonblock-cblas 2D-cannon-nonblock-cblas-prepend

clean :
	for dir in serial cannon; do \
		(cd src/$$dir; make veryclean;) \
	done;

clean-exe : clean
	rm ./x.mm_*

# DO NOT DELETE
