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
		// Eigen compatible wrapper for ScalarNode
		struct ScalarWrapper {

			ScalarWrapper() : node(ScalarNode::create_zero()) { ; }
			ScalarWrapper(int i) : node(ScalarNode::create_const(i)) { ; }
			ScalarWrapper(double d) : node(ScalarNode::create_const(d)) { ; }
			ScalarWrapper(ScalarNodePtr existing_ptr) : node(existing_ptr) { ; }

			inline ScalarWrapper operator-() const {
				return un_op<_minus_>(*this);
			}

			inline ScalarWrapper& operator+=(const ScalarWrapper& b) {
				node = bin_op<_add_>(*this, b).get();
				return *this;
			}

			inline ScalarWrapper operator+(const ScalarWrapper& b) const {
				return bin_op<_add_>(*this, b);
			}

			inline ScalarWrapper operator-(const ScalarWrapper& b) const {
				return bin_op<_sub_>(*this, b);
			}

			inline ScalarWrapper operator*(const ScalarWrapper& b) const {
				return bin_op<_mul_>(*this, b);
			}

			inline ScalarWrapper operator/(const ScalarWrapper& b) const {
				return bin_op<_div_>(*this, b);
			}

			inline bool operator==(const ScalarWrapper& b) const {
				b.get();
				return false;// return bin_op<_eq_>(*this, b); //How??????????
			}


			inline ScalarNodePtr operator*() const { //dereference
				return get();
			}

			inline ScalarNodePtr get() const {
				return node;
			}

			template<class T>
			static ScalarWrapper bin_op(const ScalarWrapper& a, const ScalarWrapper& b) {
				return ScalarWrapper(ScalarNode::create<T>(a.get(), b.get()));
			}

			template<class T>
			static ScalarWrapper un_op(const ScalarWrapper& a) {
				return ScalarWrapper(ScalarNode::create<T>(a.get()));
			}


		protected:
			ScalarNodePtr node;


		}; //ScalarWrapper

		//inline std::ostream& operator<<(std::ostream& out, const ScalarWrapper& s) {
		//
		//}

#define UN_OP(OP)											\
		inline ScalarWrapper OP(const ScalarWrapper& s) {	\
			return ScalarWrapper::un_op<_##OP##_>(s);		\
		}													\
		using std::OP;

		inline ScalarWrapper abs(const ScalarWrapper& s) {
			return ScalarWrapper::un_op<_abs_>(s); 
		}
		using std::abs;

		//https://eigen.tuxfamily.org/dox/namespaceEigen.html#a54cc34b64b4935307efc06d56cd531df
		inline ScalarWrapper abs2(const ScalarWrapper& s) {
			return s*s;
		};

		inline ScalarWrapper sqrt(const ScalarWrapper& s) {
			return ScalarWrapper::un_op<_sqrt_>(s);
		}
		using std::sqrt;

		UN_OP(log)
		UN_OP(log10)





	} //details
} //bparser

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
			RequireInitialization = 1,
			ReadCost = HugeCost,
			AddCost = HugeCost,
			MulCost = HugeCost
		};
	};
}

#endif //!INCLUDE_SCALAR_WRAPPER_HH_