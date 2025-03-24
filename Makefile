CPPFLAGS= -DNDEBUG -std=c++11 -O3
export CC=$(CXX)
INCLUDES=ext/kseq/kseq.h src/include/utils.hpp
SOURCES1=src/chainx.cpp ext/essaMEM/sparseSA.cpp  ext/essaMEM/sssort_compact.cc
SOURCES2=src/edlib_wrapper.cpp ext/edlib/edlib.cpp
SOURCES3=src/printanchors.cpp ext/essaMEM/sparseSA.cpp ext/essaMEM/sssort_compact.cc

SOURCES4=src/chainx-mininimizer.cpp

.PHONY: all test
all: chainX edlib_wrapper printanchors chainX-minimizer

chainX: $(SOURCES1) $(INCLUDES) src/include/algo.hpp src/include/parseCmdArgs.hpp
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o chainX $(SOURCES1) -lz

edlib_wrapper: $(SOURCES2) $(INCLUDES)
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o edlib_wrapper $(SOURCES2) -lz

printanchors: $(SOURCES3) $(INCLUDES)
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o printanchors $(SOURCES3) -lz

chainX-minimizer: $(SOURCES4) $(INCLUDES)
	+$(MAKE) -C ext/minimap2-2.24
	$(CXX) $(CPPFLAGS) -I ext/ -I src/include -o chainX-mininimizer $(SOURCES4) ext/minimap2-2.24/libminimap2.a -lz -lm -lpthread

test:
	# Table 1
	data/time_semiglobal/run_chainx_test.sh
	data/time_global/run_chainx_test.sh
	# Table 2
	data/correlation_semiglobal/run_chainx_test.sh MUM 20
	data/correlation_global/run_chainx_test.sh
	# Table 3
	data/correlation_semiglobal/run_chainx_test.sh MUM 20
	data/correlation_semiglobal/run_chainx_test.sh MUM 10
	data/correlation_semiglobal/run_chainx_test.sh MUM  7
	data/correlation_semiglobal/run_chainx_test.sh MEM 20
	data/correlation_semiglobal/run_chainx_test.sh MEM 10
	data/correlation_semiglobal/run_chainx_test.sh MEM  7

clean:
	+$(MAKE) -C ext/minimap2-2.24 clean
	rm -f chainX edlib_wrapper printanchors chainX-mininimizer
