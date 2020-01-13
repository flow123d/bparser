/*
 * test_tools.hh
 *
 *  Created on: Jan 10, 2020
 *      Author: jb
 */

#ifndef TEST_TEST_TOOLS_HH_
#define TEST_TEST_TOOLS_HH_


#include "assert.hh"

class AssertExc : public std::exception {
private:
    const char* expression;
    const char* file;
    int line;
    std::string message;

public:


    /// Construct an assertion failure exception
    AssertExc(const char* expression, const char* file, int line, const std::string& message = "")
        : expression(expression)
        , file(file)
        , line(line)
        , message(message)
    {}

    /// The assertion message
    virtual const char* what() const throw()
    {
        std::ostringstream outputStream;

        if (!message.empty()) {
            outputStream << message << ": ";
        }
        outputStream << "Assert: '" << expression << "'";
        outputStream << " failed in file '" << file << "' line " << line;
        std::cerr << outputStream.str();
        return outputStream.str().c_str();
    }


    ~AssertExc()
    {}
};



#ifdef NDEBUG
	#define ASSERT(EXPRESSION) ((void)0)
#else
	#define ASSERT(EXPRESSION) \
		if ( !(EXPRESSION) ) throw AssertExc(#EXPRESSION, __FILE__, __LINE__)
#endif


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
	if ( ! success) ASSERT(false);								\
}



#endif /* TEST_TEST_TOOLS_HH_ */
