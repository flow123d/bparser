/*
 * config.hh
 *
 *  Created on: Dec 29, 2019
 *      Author: jb
 */

#ifndef INCLUDE_CONFIG_HH_
#define INCLUDE_CONFIG_HH_

#ifndef DEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#endif

#ifndef NDEBUG
#ifndef DEBUG
#define DEBUG
#endif
#endif


// common declarations

typedef unsigned int uint;

// Denoting unused function parameters.
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif


#endif /* INCLUDE_CONFIG_HH_ */
