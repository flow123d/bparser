// save diagnostic state
#pragma GCC diagnostic push 

// turn off the specific warning.
#pragma GCC diagnostic ignored "-Wunused-but-set-parameter"

#undef MAX_VECTOR_SIZE
#define MAX_VECTOR_SIZE 512

#include "vectorclass.h"

// turn the warnings back on
#pragma GCC diagnostic pop