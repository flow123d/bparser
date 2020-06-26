#
#	Makefile
#
#	Created on: Dec 26, 2019
#	Author: jb
#

COMPILER         := -g++
#COMPILER        := -clang++
OPTIMIZATION_OPT :=  -O3 -mavx2
#BASE_OPTIONS     := -fmax-errors=5 #-pedantic-errors -Wall -Wextra -Werror -Wno-long-long
INCLUDES		 := -I include
BASE_OPTIONS     := $(INCLUDES) -std=c++11 -pedantic-errors -Werror=pedantic -Wall -Wextra -Werror -Wno-long-long -Wno-strict-aliasing
DBG_OPT			 := -g -DDEBUG #-fsanitize=address -static-libasan  -fno-omit-frame-pointer
#OPTIONS          := $(BASE_OPTIONS) $(OPTIMIZATION_OPT)
OPTIONS          := $(BASE_OPTIONS) $(DBG_OPT)
LINKER_OPT       := -L/usr/lib -lstdc++ -lm

ASAN_OPT         := -g -fsanitize=address -static-libasan -fno-omit-frame-pointer
MSAN_OPT         := -g -fsanitize=memory    -fno-omit-frame-pointer
LSAN_OPT         := -g -fsanitize=leak      -fno-omit-frame-pointer
USAN_OPT         := -g -fsanitize=undefined -fno-omit-frame-pointer
BUILD_SRC        := $(sort $(wildcard exprtk_*.cpp))
BUILD_LIST       := $(BUILD_SRC:%.cpp=%)


# all: $(BUILD_LIST)
# 
# $(BUILD_LIST) : %: %.cpp exprtk.hpp
# 	$(COMPILER) $(OPTIONS) -o $@ $@.cpp $(LINKER_OPT)
# 
# strip_bin :
# 	@for f in $(BUILD_LIST); do if [ -f $$f ]; then strip -s $$f; echo $$f; fi done;
# 
# valgrind :
# 	@for f in $(BUILD_LIST); do \
# 		if [ -f $$f ]; then \
# 			cmd="valgrind --leak-check=full --show-reachable=yes --track-origins=yes --log-file=$$f.log -v ./$$f"; \
# 			echo $$cmd; \
# 			$$cmd; \
# 		fi done;
# 
# pgo: exprtk_benchmark.cpp exprtk.hpp
# 	$(COMPILER) $(BASE_OPTIONS) -O3 -march=native -fprofile-generate -o exprtk_benchmark exprtk_benchmark.cpp $(LINKER_OPT)
# 	./exprtk_benchmark
# 	$(COMPILER) $(BASE_OPTIONS) -O3 -march=native -fprofile-use -o exprtk_benchmark exprtk_benchmark.cpp $(LINKER_OPT)

	
clean:
	cd build && rm -f *


	
build/grammar.o: include/grammar.hh include/grammar.cc include/grammar.impl.hh
	$(COMPILER) $(OPTIONS) -o build/grammar.o -c include/grammar.cc
	
grammar: build/grammar.o


## Proof of concept snippets.

test_design:
	rm -f build/test_design 2>/dev/null
	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_design test/test_design.cc
	$(COMPILER) $(OPTIONS) -o build/test_design test/test_design.cc
	build/test_design
	#build/test_design
	#build/test_design

#test_design:	
#	$(COMPILER) $(BASE_OPTIONS) -S -fverbose-asm -O2 -mavx2  -std=c++11 -I include  -o build/test_design.S test/test_design.cc
#	as -alhnd build/test_design.S > build/test_design.info.S

test_simd:
	rm -f build/test_simd 2>/dev/null
	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_design test/test_design.cc
	$(COMPILER) $(OPTIONS) -o build/test_simd test/test_simd.cc
	build/test_simd
	#build/test_design
	#build/test_design



## Unit tests.

test_array:
	rm -f build/test_array 2>/dev/null
	$(COMPILER) $(OPTIONS) -o build/test_array test/test_array.cc
	build/test_array
	
test_processor:
	rm -f build/test_processor 2>/dev/null
	$(COMPILER) $(OPTIONS) -o build/test_processor test/test_processor.cc
	#$(COMPILER) $(BASE_OPTIONS) -O3 -mavx2  -std=c++11 -I include  -o build/test_simd test/test_simd.cc
	build/test_processor

test_grammar: grammar
	rm -f build/test_grammar 2>/dev/null
	$(COMPILER) $(OPTIONS) -o build/test_grammar build/grammar.o test/test_grammar.cc
	build/test_grammar

test_parser: grammar
	rm -f build/test_parser 2>/dev/null
	$(COMPILER) $(OPTIONS) -o build/test_parser build/grammar.o test/test_parser.cc
	#$(COMPILER) $(BASE_OPTIONS) -O3 -mavx2  -std=c++11 -I include  -o build/test_simd test/test_simd.cc
	build/test_parser

	
test_speed: grammar
	rm -f build/test_speed 2>/dev/null
	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_speed build/grammar.o test/test_parser_speed.cc
	$(COMPILER) $(OPTIONS) -o build/test_speed build/grammar.o test/test_parser_speed.cc
	build/test_speed

tests: test_array test_processor test_parser test_speed