/*
 * test_tools.hh
 *
 *  Created on: Jan 10, 2020
 *      Author: jb
 */

#ifndef TEST_TEST_TOOLS_HH_
#define TEST_TEST_TOOLS_HH_


#include "assert.hh"



#define ASSERT_THROW(expression, msg) {				            \
	bool success = false;										\
    try {                                                       \
        expression;                                             \
    } catch (const bparser::Exception & e) {                    \
    	std::string subwhat = std::string(e.what()).substr(0, std::string(msg).size()); \
		if (msg == subwhat ) \
			success = true;									    \
		else													\
			std::cout << "Wrong exception msg: " << e.what() << "\n";	\
    } catch (const std::exception &e) {                          \
    	throw e;												\
    }                                                           \
	if ( ! success) throw;								\
}



#endif /* TEST_TEST_TOOLS_HH_ */
