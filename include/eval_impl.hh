
namespace bparser {
using namespace details;

template <typename VecType>
struct Vec {
	double *values;
	uint *subset;

	typedef VecType MyVCLVec;

	void set(double * v, uint * s) {
		values = v;
		subset = s;
	}

	inline double * value(uint i) {
		// std::cout << "self: " << this << std::endl;
		// std::cout << "v: " << values << " s: " << subset << std::endl;
		// std::cout << "i: " << i << std::endl;
		// std::cout << " si: " << subset[i] << std::endl;

		return &(values[subset[i]]);
	}


	VecType true_value() {
		return get_true_value<VecType>();
	}

	VecType false_value() {
		return get_false_value<VecType>();
	}

};


/**
 * Processor's storage.
 */
template <typename VecType>
struct Workspace {
	uint vec_n_blocks;

	// Array of vectors. Temporaries, input vectors and result vectors.
	Vec<VecType> *vector;

	uint subset_size;
	uint *const_subset;
	uint *vec_subset;

};


/**
 * Memory aligned representation of single operation.
 */
struct Operation {
	// Op code. See scalar_expr.hh: XYZNode::op_code;
	unsigned char code;
	// index of arguments in the Processors's workspace
	unsigned char arg[4];
};


template<uint NParams, class T, typename VecType>
struct EvalImpl;
//{
//	static inline void eval(Operation op, Workspace &w) {};
//};


// EvalmImpl with 1 operand
template <class T, typename VecType>
struct EvalImpl<1, T, VecType> {
	inline static void eval(Operation op, Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<1, T, double> {
	inline static void eval(Operation op, Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<1, T, VecType>::eval(Operation op, Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		VecType v0i;

		// load value into vector
		v0i.load(v0id);

		// evaluate result
		T::eval(v0i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<1, T, double>::eval(Operation op, Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		
		// evaluate result
		T::eval(*v0id);
	}
}


// EvalmImpl with 2 operands
template <class T, typename VecType>
struct EvalImpl<2, T, VecType> {
	inline static void eval(Operation op, Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<2, T, double> {
	inline static void eval(Operation op, Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<2, T, VecType>::eval(Operation op, Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		VecType v0i;
		VecType v1i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);

		// evaluate result
		T::eval(v0i, v1i);

		// store result into memory at v0id
		v0i.store(v0id); 
	}
}

template <class T>
inline void EvalImpl<2, T, double>::eval(Operation op, Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id);
	}
}


// EvalmImpl with 3 operands
template <class T, typename VecType>
struct EvalImpl<3, T, VecType> {
	inline static void eval(Operation op, Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<3, T, double> {
	inline static void eval(Operation op, Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<3, T, VecType>::eval(Operation op, Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];
	Vec<VecType> v2 = w.vector[op.arg[2]];
		// std::cout << "iv0:" << uint(op.arg[0])
		// 		<< "iv1:" << uint(op.arg[1])
		// 		<< "iv2:" << uint(op.arg[2]) << std::endl;

	for(uint i=0; i<w.subset_size; ++i) {
		// std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		VecType v0i;
		VecType v1i;
		VecType v2i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);
		v2i.load(v2id);

		// evaluate result
		T::eval(v0i, v1i, v2i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<3, T, double>::eval(Operation op, Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];
	Vec<double> v2 = w.vector[op.arg[2]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id, *v2id);
	}
}


// EvalmImpl with 4 operands
template <class T, typename VecType>
struct EvalImpl<4, T, VecType> {
	inline static void eval(Operation op, Workspace<VecType> &w);
};

template <class T>
struct EvalImpl<4, T, double> {
	inline static void eval(Operation op, Workspace<double> &w);
};

template <class T, typename VecType>
inline void EvalImpl<4, T, VecType>::eval(Operation op, Workspace<VecType> &w) {
	Vec<VecType> v0 = w.vector[op.arg[0]];
	Vec<VecType> v1 = w.vector[op.arg[1]];
	Vec<VecType> v2 = w.vector[op.arg[2]];
	Vec<VecType> v3 = w.vector[op.arg[3]];
		// std::cout << "iv0:" << uint(op.arg[0])
		// 		<< "iv1:" << uint(op.arg[1])
		// 		<< "iv2:" << uint(op.arg[2])
		// 		<< "iv3:" << uint(op.arg[3]) << std::endl;

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		double * v3id = v3.value(i);
		VecType v0i;
		VecType v1i;
		VecType v2i;
		VecType v3i;

		// load values into vectors
		v0i.load(v0id);
		v1i.load(v1id);
		v2i.load(v2id);
		v3i.load(v3id);

		// evaluate result
		T::eval(v0i, v1i, v2i, v3i);

		// store result into memory at v0id
		v0i.store(v0id);
	}
}

template <class T>
inline void EvalImpl<4, T, double>::eval(Operation op, Workspace<double> &w) {
	Vec<double> v0 = w.vector[op.arg[0]];
	Vec<double> v1 = w.vector[op.arg[1]];
	Vec<double> v2 = w.vector[op.arg[2]];
	Vec<double> v3 = w.vector[op.arg[3]];

	for(uint i=0; i<w.subset_size; ++i) {
		//std::cout << "subset: " << i << std::endl;

		double * v0id = v0.value(i);
		double * v1id = v1.value(i);
		double * v2id = v2.value(i);
		double * v3id = v3.value(i);
		
		// evaluate result
		T::eval(*v0id, *v1id, *v2id, *v3id);
	}
}

} // bparser namespace
