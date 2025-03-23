
/*
 * aligned_alloc.hh
 *
 *  Created on: Mar 23, 2025
 *      Author: LV
 */

 /*
 * Wraps the third party library function for DLL export reasons.
 */

#ifndef INCLUDE_INSTRSET_DETECT_HH
#define INCLUDE_INSTRSET_DETECT_HH

#include "config.hh"
#include "instrset.h"

EXPORT int b_instrset_detect(void);

#endif