
/*
 * aligned_alloc.hh
 *
 *  Created on: Mar 23, 2025
 *      Author: LV
 */

 /*
 * Wraps the third party library function for DLL export reasons.
 */

#ifndef INCLUDE_INSTRSET_DETECT_HH_
#define INCLUDE_INSTRSET_DETECT_HH_

#include "config.hh"
#include "instrset.h"
namespace bparser{

	EXPORT int b_instrset_detect(void);

}
#endif //!INCLUDE_INSTRSET_DETECT_HH_