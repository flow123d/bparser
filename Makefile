#
#	Makefile
#
#	Created on: Dec 26, 2019
#	Author: jb
#


	
clean:
	cd build && rm -rf *

cmake_dbg:
	cd build && cmake -DCMAKE_BUILD_TYPE=debug ..

cmake_rel:
	cd build && cmake -DCMAKE_BUILD_TYPE=release ..
	
test_array: 
	cd build && make test_array 

test_processor: 
	cd build && make test_processor 

test_grammar: 
	cd build && make test_grammar 

test_parser: 
	cd build && make test_parser
	
test_speed: 
	cd build && make test_speed
	 

#test_design:
#	rm -f build/test_design 2>/dev/null
#	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_design test/test_design.cc
#	$(COMPILER) $(OPTIONS) -o build/test_design test/test_design.cc
#	build/test_design
#	#build/test_design
#	#build/test_design
#
##test_design:	
##	$(COMPILER) $(BASE_OPTIONS) -S -fverbose-asm -O2 -mavx2  -std=c++11 -I include  -o build/test_design.S test/test_design.cc
##	as -alhnd build/test_design.S > build/test_design.info.S
#
#test_simd:
#	rm -f build/test_simd 2>/dev/null
#	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_design test/test_design.cc
#	$(COMPILER) $(OPTIONS) -o build/test_simd test/test_simd.cc
#	build/test_simd
#	#build/test_design
#	#build/test_design
#
#
#
### Unit tests.
#
#test_array:
#	rm -f build/test_array 2>/dev/null
#	$(COMPILER) $(OPTIONS) -o build/test_array test/test_array.cc
#	build/test_array
#	
#test_processor:
#	rm -f build/test_processor 2>/dev/null
#	$(COMPILER) $(OPTIONS) -o build/test_processor test/test_processor.cc
#	#$(COMPILER) $(BASE_OPTIONS) -O3 -mavx2  -std=c++11 -I include  -o build/test_simd test/test_simd.cc
#	build/test_processor
#
#test_grammar: grammar
#	rm -f build/test_grammar 2>/dev/null
#	$(COMPILER) $(OPTIONS) -o build/test_grammar build/grammar.o test/test_grammar.cc
#	build/test_grammar
#
#test_parser: grammar
#	rm -f build/test_parser 2>/dev/null
#	$(COMPILER) $(OPTIONS) -o build/test_parser build/grammar.o test/test_parser.cc
#	#$(COMPILER) $(BASE_OPTIONS) -O3 -mavx2  -std=c++11 -I include  -o build/test_simd test/test_simd.cc
#	build/test_parser
#
#	
#test_speed: grammar
#	rm -f build/test_speed 2>/dev/null
#	#$(COMPILER) $(BASE_OPTIONS) $(DBG_OPT)  -std=c++11 -I include  -o build/test_speed build/grammar.o test/test_parser_speed.cc
#	$(COMPILER) $(OPTIONS) -o build/test_speed build/grammar.o test/test_parser_speed.cc
#	build/test_speed
#
#tests: test_grammar test_array  test_processor test_parser test_speed