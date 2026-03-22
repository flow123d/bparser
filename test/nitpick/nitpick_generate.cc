#include "nitpick_include.hh"
#include <fstream>
#include <filesystem>
#define NITPICK_IDE_IGNORE

#include NITPICK_DEF_FILE

using namespace bparser;

int main() {

	const uint max_vec_size = 1;
	const size_t buffer_size = sizeof(double) * 2000;

	void* buffer = ::operator new(buffer_size);
	

	//double buffer[buffer_size];
	ExprCase exprcase(max_vec_size,(void*)buffer, buffer_size);

	//def is in NITPICK_DEF_FILE
	def(exprcase);

	if (exprcase.get_expression().empty()) {
		Throw() << "No expression was set!";
	}

	exprcase.allocate(1);
	Parser& p = exprcase.get_parser();

	// Compile the expression into internal processor.
	p.compile();
	//p.compile(exprcase.get_patch_arena()); //Add arena from ExprCase

	ExpressionDAG dag(p.result_array().elements());
	DagPrinter(dag).print_in_dot2(p.get_raw_symbols());

	std::ofstream file(NITPICK_GEN_FILE);
	file << DagPrinter(dag).print_in_cxx(exprcase.get_inv_map());
	file.close();
	std::cout << "File " << std::filesystem::current_path() << " " << NITPICK_GEN_FILE << " created from " << NITPICK_DEF_FILE;

	::operator delete(buffer);

} //main