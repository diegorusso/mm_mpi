default : all

serial : 
	if test -d mm-serial ; then \
	( cd mm-serial ; make ; cp x.mm ../x.mm_serial ; ) ; fi

2D-cannon : 
	if test -d mm-2D-cannon ; then \
	( cd mm-2D-cannon ; make ; cp x.mm ../x.mm_2D_cannon ; ) ; fi

2D-cannon-cblas : 
	if test -d mm-2D-cannon ; then \
	( cd mm_2D-cannon ; make clean; make -f Makefile.cblas ; cp x.mm ../x.mm_2D_cannon_cblas ; ) ; fi

2D-cannon-nonblock : 
	if test -d mm-2D-cannon-nonblock ; then \
	( cd mm-2D-cannon-nonblock ; make ; cp x.mm ../x.mm_2D_cannon_nonblock ; ) ; fi

2D-cannon-nonblock-cblas : 
	if test -d mm-2D-cannon-nonblock ; then \
	( cd mm-2D-cannon-nonblock ; make clean ; make -f Makefile.cblas ; cp x.mm ../x.mm_2D_cannon_cblas_nonblock ; ) ; fi

2D-cannon-nonblock-cblas-prepend : 
	if test -d mm-2D-cannon-nonblock ; then \
	( cd mm-2D-cannon-nonblock ; make clean ; make -f Makefile.cblas-prepend ; cp x.mm ../x.mm_2D_cannon_cblas_nonblock_prepend ; ) ; fi

all : serial 2D-cannon 2D-cannon-nonblock 2D-cannon-cblas 2D-cannon-nonblock-cblas 2D-cannon-nonblock-cblas-prepend

clean :
	for dir in \
		2D-cannon 2D-cannon-nonblock 2D-cannon-cblas 2D-cannon-nonblock-cblas 2D-cannon-nonblock-cblas-prepend \
	; do \
	    if test -d mm-$$dir ; then \
		( cd mm-$$dir ; \
		make veryclean ; ) \
	    fi \
	done ;  \
	if test -d mm ; then \
		( cd mm ; \
		make veryclean ; ) \
	fi

clean-exe : clean
	rm ./x.mm_*

# DO NOT DELETE
