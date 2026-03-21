#include "nitpick_include.hh"
#include "testcase.hh"
//#include "scalar_node.hh"


using namespace bparser;


int main() {
	

	//For ExprCase arena
	const size_t buffer_size = sizeof(double) * TestCase::N * 32;

	std::cout << "Buffer size: " << buffer_size << std::endl;
	{
		int simd_size = bparser::get_simd_size();
		std::cout << "SIMD size: " << simd_size << std::endl;
	}

	void* buffer = ::operator new(buffer_size);
		
	std::vector<uint> vec_sizes({ 16U, 64U, 256U, 1024U, TestCase::N });


	for (uint vec_size : vec_sizes) {

		EfficientTestCaseForceVariable(vec_size, buffer, buffer_size).run();
		ShiftingTestCaseForceVariable(vec_size, buffer, buffer_size).run();
		//ShiftingTestCaseForceCopy(vec_size, buffer, buffer_size).run();
		
	} //for


	//std::cout << "Result: \n";

	// Result in the 'vres' value space.
	//std::cout << print_vec(vres, vres_size);
	

	::operator delete(buffer);

} //main