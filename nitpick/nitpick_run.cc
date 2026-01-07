#include "nitpick_include.hh"

int main() {
// only get us the variables
#define NITPICK_IDE_IGNORE
#include "nitpick_common.cc" 

	//include all the variables used to generate the file
#include NITPICK_DEF_FILE
	//include the generated file
#include NITPICK_GEN_FILE

	//se is in the NITPICK_GEN_FILE
	//max_vec_size is in NITPICK_DEF_FILE
	ProcessorBase* processor = ProcessorBase::create_processor(se, max_vec_size, bparser::get_simd_size(), nullptr);

	std::vector<uint> subset = { 0, 1 }; //ctverice doubluu //TODO: move or expand

	processor->set_subset(subset);
	processor->run();

	//TODO:
	std::cout << "Result: \n";

	// Result in the 'vres' value space.
	std::cout << print_vec(vres, vres_size);
}