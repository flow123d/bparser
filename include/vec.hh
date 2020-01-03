#ifndef Vec_HH
#define Vec_HH


#include "assert.hh"
#include "numeric.hh"


/**
 * Vec - is a basic data type of the bparser.
 *
 * Workspace - kind of custom allocator of Vec contents and constants.
 * workspace is always aligned to the size of the vector operations.
 *
 *
 * Currently only double is supported.
 */

namespace bparser {
struct Vec;


/**
 *
 */
class Workspace {
public:
	// Size of the single vector operation. E.g. 4 doubles for AVX2.
	static const uint simd_block_size = 4;

	Workspace()
	: Vec_size_(0)
	{
		workspace_size_ = 64;
		workspace_ = new double[workspace_size_];
		clear();
	}

	~Workspace() {
		delete [] workspace_;
	}


	static void set_subset(std::initializer_list<uint> subset, uint Vec_size) {
		Workspace::instance().set_subset_(std::vector<uint>(subset), Vec_size);
	}

	// Set new subset structure and size of the full Vec.
	static void set_subset(const std::vector<uint> &subset, uint Vec_size) {
		Workspace::instance().set_subset_(subset, Vec_size);
	}

	static void set_workspace(uint n_doubles) {
		Workspace::instance().set_workspace_(n_doubles);
	}

	// Release all temporary slots.
	static void clear() {
		Workspace::instance().next_slot_ = Workspace::instance().workspace_;
	}

	static uint size() {
		return Workspace::instance().subset_.size();
	}

	static double * get_slot()  {
		return Workspace::instance().get_slot_();
	}


protected:
	static inline Workspace &instance() {
		static Workspace w;
		return w;
	}

	inline void set_subset_(const std::vector<uint> &subset, uint Vec_size) {
		clear();
		Vec_size_ = Vec_size;
		subset_ = subset; // TODO: avoid copy
		const_subset_.reserve(size());
		flat_subset_.reserve(size());
		for(uint i=0; i<size(); ++i) {
			const_subset_[i] = 0;
			flat_subset_[i] = i;
		}
	}

	void set_workspace_(uint n_doubles) {
		delete [] workspace_;
		workspace_size_ = n_doubles;
		workspace_ = new double[n_doubles];
	}

	double * get_slot_()  {
		ASSERT(next_slot_ < workspace_ + workspace_size_);
		double *ptr = next_slot_;
		next_slot_ += size();
		return ptr;
	}


	static uint * subset() {
		return &(Workspace::instance().subset_[0]);
	}

	static uint * flat_subset() {
		return &(Workspace::instance().flat_subset_[0]);
	}

	static uint * const_subset() {
		return &(Workspace::instance().const_subset_[0]);
	}


	std::vector<uint > const_subset_;
	std::vector<uint > flat_subset_;
	std::vector<uint > subset_;
	uint Vec_size_;

	double *workspace_;
	uint workspace_size_;
	double *next_slot_;
	friend struct exprtk::Vec;
};


struct Iter {
	Iter() : i(0), j(0) {}
	void inc() {
		++j;
		if (j == Workspace::simd_block_size) {
			j=0;
			++i;
		}
	}
	uint i;
	uint j;
};

#define AFOR(it) \
	for(Iter it; it.i < size_; ++it.i)  \
		for(it.j=0; it.j < Workspace::simd_block_size; ++it.j)

struct Block {
	uint start;
	uint size;
};
typedef std::vector<Block> Subset;

struct Vec {
	typedef Workspace workspace_t;

	static Vec epsilon() {
		return Vec(0.0000000001);
	}

	Vec(double *ptr)
	: ptr_(ptr),
	  subset_(workspace_t::subset()),
	  size_(workspace_t::size()),
	  const_value_(0)
	{}

	Vec()
	: ptr_(workspace_t::get_slot()),
	  subset_(workspace_t::flat_subset()),
	  size_(workspace_t::size()),
	  const_value_(0)
	{}

	Vec(double x)
	: ptr_(&const_value_),
	  subset_(workspace_t::const_subset()),
	  size_(workspace_t::size()),
	  const_value_(x)
	{}

	Vec(const std::initializer_list<double> &list) {
		Iter it;
		for(auto item : list) {
			at(it) = item;
			it.inc();
		}
	}

	Vec(const Vec &other)
	: Vec()
	{ *this = other; }

	const Vec & operator=(const Vec &other) {
		ptr_ = other.ptr_;
		subset_ = other.subset_;
		size_ = other.size_;
		const_value_ = other.const_value_;
		return *this;
	}

	double & at(const Iter &it) {
		return *(ptr_ + *(subset_ + it.i) + it.j);
	}

	const double & at(const Iter &it) const {
		return *(ptr_ + *(subset_ + it.i) + it.j);
	}

	double & at(uint i) {
		return *(ptr_ + *(subset_ + i/workspace_t::simd_block_size) + i%workspace_t::simd_block_size);
	}

	uint size() const {
		return size_;
	}

	void reset_data(double *ptr) {
		ptr_ = ptr;
	}

	double * data() const {
		return ptr_;
	}

	void print(std::ostream &stream) {
		stream << "[ ";
		AFOR(it) {
			stream << at(it) << ", ";
		}
		stream << "]";
	}

private:
	double *ptr_;
	uint *subset_;
	uint size_;
};


Vec operator+(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) + b.at(it);
	return z;
}

Vec operator-(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) - b.at(it);
	return z;
}

Vec operator*(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) * b.at(it);
	return z;
}

Vec operator/(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) / b.at(it);
	return z;
}

Vec operator<(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) < b.at(it);
	return z;
}

Vec operator<=(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) <= b.at(it);
	return z;
}

Vec operator>(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) > b.at(it);
	return z;
}

Vec operator>=(const Vec & a, const Vec & b) {
	Vec z; uint size_ = a.size();
	AFOR(it) z.at(it) = a.at(it) >= b.at(it);
	return z;
}

std::ostream& operator<<(std::ostream& stream, Vec arr) {
    arr.print(stream);
	return stream;
}





} // exprtk


double dinf = std::numeric_limits<double>::infinity();
double dnan = std::numeric_limits<double>::quiet_NaN();
namespace std {
  template <>
  class numeric_limits<exprtk::Vec> {
  public:
    static exprtk::Vec infinity() {return exprtk::Vec(dinf);}
    static exprtk::Vec quiet_NaN() {return exprtk::Vec(dnan);}
  };
}





namespace exprtk {
namespace details {

namespace numeric {
namespace details {




template <typename T> inline bool is_nan_impl(const T v, Vec_type_tag) {
  return std::not_equal_to<T>()(v, v);
}

template <typename T> inline int to_int32_impl(const T v, Vec_type_tag) {
  return static_cast<int>(v);
}

template <typename T>
inline long long int to_int64_impl(const T v, Vec_type_tag) {
  return static_cast<long long int>(v);
}

//template <typename T> inline bool is_true_impl(const T v) {
//  return std::not_equal_to<T>()(T(0), v);
//}
//
//template <typename T> inline bool is_false_impl(const T v) {
//  return std::equal_to<T>()(T(0), v);
//}

template <typename T> inline T abs_impl(const T v, Vec_type_tag) {
  return ((v < T(0)) ? -v : v);
}

template <typename T> inline T min_impl(const T v0, const T v1, Vec_type_tag) {
  return std::min<T>(v0, v1);
}

template <typename T> inline T max_impl(const T v0, const T v1, Vec_type_tag) {
  return std::max<T>(v0, v1);
}

template <typename T> inline T equal_impl(const T v0, const T v1, Vec_type_tag) {
  return (abs_impl(v0 - v1, Vec_type_tag()) <=
          (std::max(T(1), std::max(abs_impl(v0, Vec_type_tag()),
                                   abs_impl(v1, Vec_type_tag()))) *
           Vec::epsilon()))
             ? T(1)
             : T(0);
}


template <typename T> inline T expm1_impl(const T v, Vec_type_tag) {
  // return std::expm1<T>(v);
  if (abs_impl(v, Vec_type_tag()) < T(0.00001))
    return v + (T(0.5) * v * v);
  else
    return std::exp(v) - T(1);
}


template <typename T>
inline T nequal_impl(const T v0, const T v1, Vec_type_tag) {
  typedef Vec_type_tag rtg;
  return (abs_impl(v0 - v1, rtg()) >
          (std::max(T(1), std::max(abs_impl(v0, rtg()), abs_impl(v1, rtg()))) *
           Vec::epsilon()))
             ? T(1)
             : T(0);
}

template <typename T>
inline T modulus_impl(const T v0, const T v1, Vec_type_tag) {
  return std::fmod(v0, v1);
}


template <typename T> inline T pow_impl(const T v0, const T v1, Vec_type_tag) {
  return std::pow(v0, v1);
}


template <typename T>
inline T logn_impl(const T v0, const T v1, Vec_type_tag) {
  return std::log(v0) / std::log(v1);
}


template <typename T> inline T log1p_impl(const T v, Vec_type_tag) {
  if (v > T(-1)) {
    if (abs_impl(v, Vec_type_tag()) > T(0.0001)) {
      return std::log(T(1) + v);
    } else
      return (T(-0.5) * v + T(1)) * v;
  } else
    return std::numeric_limits<T>::quiet_NaN();
}


template <typename T>
inline T root_impl(const T v0, const T v1, Vec_type_tag) {
  if (v1 < T(0))
    return std::numeric_limits<T>::quiet_NaN();

  const std::size_t n = static_cast<std::size_t>(v1);

  if ((v0 < T(0)) && (0 == (n % 2)))
    return std::numeric_limits<T>::quiet_NaN();

  return std::pow(v0, T(1) / n);
}


template <typename T> inline T round_impl(const T v, Vec_type_tag) {
  return ((v < T(0)) ? std::ceil(v - T(0.5)) : std::floor(v + T(0.5)));
}

//template <typename T>
//inline T roundn_impl(const T v0, const T v1, Vec_type_tag) {
//  const int index =
//      std::max<int>(0, std::min<int>(pow10_size - 1, (int)std::floor(v1)));
//  const T p10 = T(pow10[index]);
//
//  if (v0 < T(0))
//    return T(std::ceil((v0 * p10) - T(0.5)) / p10);
//  else
//    return T(std::floor((v0 * p10) + T(0.5)) / p10);
//}


template <typename T>
inline T hypot_impl(const T v0, const T v1, Vec_type_tag) {
  return std::sqrt((v0 * v0) + (v1 * v1));
}


template <typename T>
inline T atan2_impl(const T v0, const T v1, Vec_type_tag) {
  return std::atan2(v0, v1);
}


template <typename T> inline T shr_impl(const T v0, const T v1, Vec_type_tag) {
  return v0 * (T(1) / std::pow(T(2), static_cast<T>(static_cast<int>(v1))));
}


template <typename T> inline T shl_impl(const T v0, const T v1, Vec_type_tag) {
  return v0 * std::pow(T(2), static_cast<T>(static_cast<int>(v1)));
}


template <typename T> inline T sgn_impl(const T v, Vec_type_tag) {
  if (v > T(0))
    return T(+1);
  else if (v < T(0))
    return T(-1);
  else
    return T(0);
}

template <typename T> inline T and_impl(const T v0, const T v1, Vec_type_tag) {
  return (is_true_impl(v0) && is_true_impl(v1)) ? T(1) : T(0);
}


template <typename T>
inline T nand_impl(const T v0, const T v1, Vec_type_tag) {
  return (is_false_impl(v0) || is_false_impl(v1)) ? T(1) : T(0);
}


template <typename T> inline T or_impl(const T v0, const T v1, Vec_type_tag) {
  return (is_true_impl(v0) || is_true_impl(v1)) ? T(1) : T(0);
}


template <typename T> inline T nor_impl(const T v0, const T v1, Vec_type_tag) {
  return (is_false_impl(v0) && is_false_impl(v1)) ? T(1) : T(0);
}


template <typename T> inline T xor_impl(const T v0, const T v1, Vec_type_tag) {
  return (is_false_impl(v0) != is_false_impl(v1)) ? T(1) : T(0);
}


template <typename T>
inline T xnor_impl(const T v0, const T v1, Vec_type_tag) {
  const bool v0_true = is_true_impl(v0);
  const bool v1_true = is_true_impl(v1);

  if ((v0_true && v1_true) || (!v0_true && !v1_true))
    return T(1);
  else
    return T(0);
}



template <typename T> inline T ncdf_impl(T v, Vec_type_tag) {
  T cnd = T(0.5) * (T(1) + erf_impl(abs_impl(v, Vec_type_tag()) /
                                        T(numeric::constant::sqrt2),
                                    Vec_type_tag()));
  return (v < T(0)) ? (T(1) - cnd) : cnd;
}

template <typename T> inline T sinc_impl(T v, Vec_type_tag) {
  if (std::abs(v) >= std::numeric_limits<T>::epsilon())
    return (std::sin(v) / v);
  else
    return T(1);
}


template <typename T> inline T acos_impl(const T v, Vec_type_tag) {
  return std::acos(v);
}
template <typename T> inline T acosh_impl(const T v, Vec_type_tag) {
  return std::log(v + std::sqrt((v * v) - T(1)));
}
template <typename T> inline T asin_impl(const T v, Vec_type_tag) {
  return std::asin(v);
}
template <typename T> inline T asinh_impl(const T v, Vec_type_tag) {
  return std::log(v + std::sqrt((v * v) + T(1)));
}
template <typename T> inline T atan_impl(const T v, Vec_type_tag) {
  return std::atan(v);
}
template <typename T> inline T atanh_impl(const T v, Vec_type_tag) {
  return (std::log(T(1) + v) - std::log(T(1) - v)) / T(2);
}
template <typename T> inline T ceil_impl(const T v, Vec_type_tag) {
  return std::ceil(v);
}
template <typename T> inline T cos_impl(const T v, Vec_type_tag) {
  return std::cos(v);
}
template <typename T> inline T cosh_impl(const T v, Vec_type_tag) {
  return std::cosh(v);
}
template <typename T> inline T exp_impl(const T v, Vec_type_tag) {
  return std::exp(v);
}
template <typename T> inline T floor_impl(const T v, Vec_type_tag) {
  return std::floor(v);
}
template <typename T> inline T log_impl(const T v, Vec_type_tag) {
  return std::log(v);
}
template <typename T> inline T log10_impl(const T v, Vec_type_tag) {
  return std::log10(v);
}
template <typename T> inline T log2_impl(const T v, Vec_type_tag) {
  return std::log(v) / T(numeric::constant::log2);
}
template <typename T> inline T neg_impl(const T v, Vec_type_tag) { return -v; }
template <typename T> inline T pos_impl(const T v, Vec_type_tag) { return +v; }
template <typename T> inline T sin_impl(const T v, Vec_type_tag) {
  return std::sin(v);
}
template <typename T> inline T sinh_impl(const T v, Vec_type_tag) {
  return std::sinh(v);
}
template <typename T> inline T sqrt_impl(const T v, Vec_type_tag) {
  return std::sqrt(v);
}
template <typename T> inline T tan_impl(const T v, Vec_type_tag) {
  return std::tan(v);
}
template <typename T> inline T tanh_impl(const T v, Vec_type_tag) {
  return std::tanh(v);
}
template <typename T> inline T cot_impl(const T v, Vec_type_tag) {
  return T(1) / std::tan(v);
}
template <typename T> inline T sec_impl(const T v, Vec_type_tag) {
  return T(1) / std::cos(v);
}
template <typename T> inline T csc_impl(const T v, Vec_type_tag) {
  return T(1) / std::sin(v);
}
template <typename T> inline T r2d_impl(const T v, Vec_type_tag) {
  return (v * T(numeric::constant::_180_pi));
}
template <typename T> inline T d2r_impl(const T v, Vec_type_tag) {
  return (v * T(numeric::constant::pi_180));
}
template <typename T> inline T d2g_impl(const T v, Vec_type_tag) {
  return (v * T(20.0 / 9.0));
}
template <typename T> inline T g2d_impl(const T v, Vec_type_tag) {
  return (v * T(9.0 / 20.0));
}
template <typename T> inline T notl_impl(const T v, Vec_type_tag) {
  return (std::not_equal_to<T>()(T(0), v) ? T(0) : T(1));
}
template <typename T> inline T frac_impl(const T v, Vec_type_tag) {
  return (v - static_cast<long long>(v));
}
template <typename T> inline T trunc_impl(const T v, Vec_type_tag) {
  return T(static_cast<long long>(v));
}

template <typename T> inline T const_pi_impl(Vec_type_tag) {
  return T(numeric::constant::pi);
}
template <typename T> inline T const_e_impl(Vec_type_tag) {
  return T(numeric::constant::e);
}

template <typename T> inline bool is_integer_impl(const T &v, Vec_type_tag) {
  return std::equal_to<T>()(T(0), std::fmod(v, T(1)));
}

} // namespace details
} // namespace numeric



} // namespace details
} // namespace exprtk


#endif
