#include "nitpick_include.hh"
#include <fstream>
#include <filesystem>

int main() {

#define NITPICK_IDE_IGNORE
#include "nitpick_common.cc"

#include NITPICK_DEF_FILE

	//p is in NITPICK_DEF_FILE
	// Compile the expression into internal processor.
	p.compile();

	ExpressionDAG dag(p.result_array().elements());
	//dag.print_in_dot2();

	std::ofstream file(NITPICK_GEN_FILE);
	file << DagPrinter(dag).print_in_cxx(inv_map);
	file.close();
	std::cout << "File " << std::filesystem::current_path() << " " << NITPICK_GEN_FILE << " created from " << NITPICK_DEF_FILE;
}