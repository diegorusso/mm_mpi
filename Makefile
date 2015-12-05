default: all

## Serial version
serial:
	cd src/serial && make && cp x.mm ../../x.mm_serial;

## Cannon - blocking version
2D-cannon:
	cd src/cannon && make -f Makefile && cp x.mm ../../x.mm_2D_cannon;

# Optimizations
2D-cannon-openmp-outer:
	cd src/cannon && make -f Makefile.openmp.outer && cp x.mm ../../x.mm_2D_cannon_openmp_outer;

2D-cannon-openmp-middle:
	cd src/cannon && make -f Makefile.openmp.middle && cp x.mm ../../x.mm_2D_cannon_openmp_middle;

2D-cannon-openmp-inner:
	cd src/cannon && make -f Makefile.openmp.inner && cp x.mm ../../x.mm_2D_cannon_openmp_inner;

2D-cannon-openmp-nested:
	cd src/cannon && make -f Makefile.openmp.nested && cp x.mm ../../x.mm_2D_cannon_openmp_nested;

2D-cannon-cblas:
	cd src/cannon && make -f Makefile.cblas && cp x.mm ../../x.mm_2D_cannon_cblas;

## Cannon - non blocking version
2D-cannon-nonblock:
	cd src/cannon && make -f Makefile.nonblock && cp x.mm ../../x.mm_2D_cannon_nonblock;

# Optimizations
2D-cannon-nonblock-openmp-outer:
	cd src/cannon && make -f Makefile.nonblock.openmp.outer && cp x.mm ../../x.mm_2D_cannon_nonblock_openmp_outer;

2D-cannon-nonblock-openmp-middle:
	cd src/cannon && make -f Makefile.nonblock.openmp.middle && cp x.mm ../../x.mm_2D_cannon_nonblock_openmp_middle;

2D-cannon-nonblock-openmp-inner:
	cd src/cannon && make -f Makefile.nonblock.openmp.inner && cp x.mm ../../x.mm_2D_cannon_nonblock_openmp_inner;

2D-cannon-nonblock-openmp-nested:
	cd src/cannon && make -f Makefile.nonblock.openmp.nested && cp x.mm ../../x.mm_2D_cannon_nonblock_openmp_nested;

2D-cannon-nonblock-cblas:
	cd src/cannon && make -f Makefile.nonblock.cblas && cp x.mm ../../x.mm_2D_cannon_nonblock_cblas;

all: serial \
	 2D-cannon 2D-cannon-openmp-outer 2D-cannon-openmp-middle 2D-cannon-openmp-inner 2D-cannon-openmp-nested 2D-cannon-cblas \
	 2D-cannon-nonblock 2D-cannon-nonblock-openmp-outer 2D-cannon-nonblock-openmp-middle 2D-cannon-nonblock-openmp-inner \
	 	2D-cannon-nonblock-openmp-nested 2D-cannon-nonblock-cblas

clean:
	cd src/serial && make clean;
	cd src/cannon && make -f Makefile clean;
	cd src/cannon && make -f Makefile.nonblock clean;
	rm ./x.mm_*
