#
#	Makefile
#
#	Created on: Dec 26, 2019
#	Author: jb
#

COMPILER         := -c++
#COMPILER        := -clang++
OPTIMIZATION_OPT := -O1
BASE_OPTIONS     := -fmax-errors=5 #-pedantic-errors -Wall -Wextra -Werror -Wno-long-long
OPTIONS          := $(BASE_OPTIONS) $(OPTIMIZATION_OPT)
LINKER_OPT       := -L/usr/lib -lstdc++ -lm
DBG_OPT			 := -g -DDEBUG -fsanitize=address   -fno-omit-frame-pointer

ASAN_OPT         := -g -fsanitize=address   -fno-omit-frame-pointer
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
	rm -f core.* *~ *.o *.bak *stackdump gmon.out *.gcda *.gcno *.gcnor *.gch
	
test_grammar:
	$(COMPILER) $(BASE_OPTIONS) $(ASAN_OPT)  -std=c++11 -I include  -o build/test_grammar test/test_grammar.cc
	build/test_grammar
	
#test_design:
#	rm -f build/test_design 2>/dev/null
#	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_design test/test_design.cc
#	$(COMPILER) $(BASE_OPTIONS) -O2 -mavx2  -std=c++11 -I include  -o build/test_design test/test_design.cc
#	build/test_design
#	#build/test_design
#	#build/test_design

test_design:	
	$(COMPILER) $(BASE_OPTIONS) -S -fverbose-asm -O2 -mavx2  -std=c++11 -I include  -o build/test_design.S test/test_design.cc
	as -alhnd build/test_design.S > build/test_design.info.S
	

	
# test_loop:
# 	#$(COMPILER) $(BASE_OPTIONS) -g -std=c++11 -I include test_array_loop.cc
# 	#$(COMPILER) $(BASE_OPTIONS) -O2 -std=c++11 -I include test_array_loop.cc
# 	$(COMPILER) $(BASE_OPTIONS) -O3 -std=c++11 -mavx2 -ffast-math -ftree-vectorizer-verbose=10 -I include test_array_loop.cc
# 	#$(COMPILER) $(BASE_OPTIONS) -O3 -std=c++11 -mavx2  -ftree-vectorizer-verbose=10 -I include test_array_loop.cc
# 	#$(COMPILER) $(BASE_OPTIONS) -O3 -march=native -std=c++11 -I include test_array_loop.cc
# 	#$(COMPILER) $(BASE_OPTIONS) -O3 -mavx2 -std=c++11 -I include test_array_loop.cc
# 	./a.out
# 	./a.out
# 	./a.out
