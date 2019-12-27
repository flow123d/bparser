#ifndef ASSERT_HH
#define ASSERT_HH

#include <exception>
#include <string>
#include <sstream>
#include <iostream>


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


#endif
