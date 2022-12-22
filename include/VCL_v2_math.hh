//test if gnuc
#ifdef __GNUC__

// save diagnostic state
#pragma GCC diagnostic push 

// turn off the specific warning.
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

// end if
#endif

// #include "vectorclass.h"

#include "vectormath_exp.h"
#include "vectormath_hyp.h"
#include "vectormath_trig.h"

//test if gnuc
#ifdef __GNUC__

// turn the warnings back on
#pragma GCC diagnostic pop

// end if
#endif
