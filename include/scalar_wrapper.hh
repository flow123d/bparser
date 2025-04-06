/*
 * scalar_wrapper.hh
 *
 *  Created on: Apr 6, 2025
 *      Author: LV
 */

//https://eigen.tuxfamily.org/dox/TopicCustomizing_CustomScalar.html

#ifndef INCLUDE_SCALAR_WRAPPER_HH_
#define INCLUDE_SCALAR_WRAPPER_HH_

#include "scalar_node.hh"
#include <Eigen/Core>

namespace bparser {
	namespace details {
		struct ScalarWrapper {

			ScalarWrapper(ScalarNodePtr existing_ptr) : node(existing_ptr) { ; }


			inline ScalarWrapper operator+(const ScalarWrapper& b) const {
				return ScalarWrapper(ScalarNode::create<_add_>(node, b.node));
			}

			inline ScalarNodePtr operator*() const {
				return get();
			}

			inline ScalarNodePtr get() const {
				return node;
			}



		protected:
			ScalarNodePtr node;
		};


	}
}
//https://eigen.tuxfamily.org/dox/structEigen_1_1NumTraits.html
namespace Eigen {
	template<> struct NumTraits<bparser::details::ScalarWrapper>
	: NumTraits<double>
	{
		typedef bparser::details::ScalarWrapper Real;
		typedef bparser::details::ScalarWrapper NonInteger;
		typedef bparser::details::ScalarWrapper Nested;

		enum {
			IsComplex = 0,
			IsInteger = 0,
			IsSigned = 1,
			RequireInitialization = 0,
			ReadCost = HugeCost,
			AddCost = HugeCost,
			MulCost = HugeCost
		};
	};
}

#endif //!INCLUDE_SCALAR_WRAPPER_HH_