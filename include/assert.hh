#ifndef ASSERT_HH
#define ASSERT_HH

#include <exception>
#include <string>
#include <sstream>
#include <iostream>

namespace bparser {

class AssertExc : public std::exception {
private:
    std::string expression;
    std::string file;
    int line;
    std::string message;

public:


    /// Construct an assertion failure exception
    AssertExc(const char* expression, const char* file, int line, const std::string& message = "")
        : expression(expression)
        , file(file)
        , line(line)
        , message(message)
    {

     }

    /// The assertion message
    virtual const char* what() const throw()
    {
    	static std::string out_message(1024,' ');

        std::ostringstream outputStream;

        if (!message.empty()) {
            outputStream << message << ": ";
        }
        outputStream << "Assert: '" << expression << "'";
        outputStream << " failed in file '" << file << "' line " << line;
        std::cerr << outputStream.str();
        out_message =  outputStream.str();
        return out_message.c_str();
    }


    ~AssertExc()
    {}
};



#ifdef NDEBUG
	#define BP_ASSERT(EXPRESSION) ((void)0)
#else
	#define BP_ASSERT(EXPRESSION) \
		if ( !(EXPRESSION) ) throw ::bparser::AssertExc(#EXPRESSION, __FILE__, __LINE__)
#endif


class Exception: public std::exception
{
public:
  mutable std::ostringstream out;
  Exception(const char *file, int line)
  {

	  out << "Bparser, " << file << ":" << line << std::endl;
	  out << "Error: ";
  }

  Exception(const Exception &other)
  : out(other.out.str())
  {}

  ~Exception() throw() {}

  virtual const char* what() const throw()
  {
     // have preallocated some space for error message we want to return
     // Is there any difference, if we move this into ExceptionBase ??
     static std::string message(1024,' ');

     // Be sure that this function do not throw.
     try {
         message = out.str();
         return message.c_str();

     } catch (std::exception &exc) {
         std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Message is " << std::endl
                 << exc.what() << std::endl << "*** Aborting! ***" << std::endl;
         std::abort();
     } catch (...) {
         std::cerr << "*** Exception encountered in exception handling routines ***" << std::endl << "*** Aborting! ***"
                 << std::endl;
         std::abort();
     }
     return 0;


  }
};

template <class T>
const Exception &operator<<(const Exception &ex, const T &v) {
	ex.out << v;
	return ex;
}

#define Throw(arg) throw Exception(__FILE__, __LINE__)

}

#endif
