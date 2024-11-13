#INCLUDES = syspara.h

bin_PROGRAMS = pvc2024 

pvc2024_SOURCES = ./src/UHNpvc.c ./src/exptable.c ./src/mem.c ./src/val_constants.c ./src/system.c ./src/currents.c\
	./src/eular.c ./src/runge.c ./src/win_out.c ./src/eventloop.c ./src/input_para.c ./src/data_out.c ./src/syspara.h

#test_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -openmp -Winline -inline-level=2 -ip -opt-mem-layout-trans=3 -pg -par-report2
#test_CFLAGS = -O2 -Wall -axCORE-AVX2 -mkl -parallel -Winline -inline-level=2 -ip -qopt-mem-layout-trans=3 -pg -qopt-report-phase=all
pvc2024_CFLAGS = -fcommon -fiopenmp -fopenmp-targets=spir64 -flto -Wall -O2 -xCORE-AVX2 -qopt-zmm-usage=high -qmkl=parallel -Winline -ipo -pg -qopt-report-stdout -parallel-source-info=2

#test_CFLAGS = -O2 -axsse4.2 -parallel -qopt-report-phase=all -qopt-report=2

#LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
pvc2024_LDADD = /home/tsumoto/lib/icx2024/xhplot.o

