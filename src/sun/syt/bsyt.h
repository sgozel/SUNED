// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_BSYT_H
#define SUN_BSYT_H

#include <iostream>
#include <string>
#include <climits>
#include <stdexcept>

#include "../../common/datatypes.h"

namespace sun {

template<class type_t, unsigned int extent>
class bSYT
{

public:
	
    bSYT() {v_=0x0ULL;};
    explicit bSYT(const type_t v) {v_ = v;};

    type_t value() const {return v_;};
    
    using value_type = type_t;
    
    static constexpr unsigned int bitextent() {return extent;}
    
    static void check(const unsigned int N, const unsigned int n) {
		constexpr unsigned int actual_maxN = (1U << extent);
		const unsigned int actual_nbits = CHAR_BIT * sizeof(type_t);
		unsigned int required_extent = 0;
		{
			unsigned int val = N - 1;
			while (val > 0) {
				val >>= 1;
				required_extent += 1;
			}
		}
		const unsigned int required_nbits = n * required_extent;
		
		if (N > actual_maxN) {
			throw std::runtime_error(std::string("Not enough bits per box for SU(") + std::to_string(N) + ") " + 
									 "[ actual: " + std::to_string(extent) + " | needed: " +
									 std::to_string(required_extent) + " ]");
		}
		if (required_nbits > actual_nbits) {
			throw std::runtime_error(std::string("Not enough bits in the container ") + 
									 "[ actual: " + std::to_string(actual_nbits) + 
									 " | needed: " + std::to_string(required_nbits) + " ]");
		}
		if (required_extent < extent) {
			std::cerr << "WARNING: bit extent is suboptimal for SU(" << std::to_string(N) << ") "
					  << "[ actual: " << std::to_string(extent)
					  << " | needed: " << std::to_string(required_extent) << " ]" << std::endl;
		}
		if (required_nbits <= actual_nbits/2) {
			std::cerr << "WARNING: container is suboptimal for SU(" << std::to_string(N) << ") with "
					  << std::to_string(n) << " boxes [ actual: " << std::to_string(actual_nbits) 
					  << " | needed " << std::to_string(required_nbits) << " ]" << std::endl;
		}
	}

	template<typename T>
    int get(const T i) const {
        int r = (v_>>(i*extent))&gmask0();
        return r;
    }
	
	template <typename T, typename U>
	void set(const T i, const U value) {
        const T pos = i*extent;
        const type_t mask = ~(((1ULL<<extent)-1)<<pos);
        v_ = (v_&mask)|(static_cast<type_t>(value)<<pos);
    }
    
    template <typename T>
	void exchange(const T i, const T j) {
		const int posi = i*extent;
		const int posj = j*extent;
		const UINT32 mask = gmask0();
		const UINT32 set1 = (v_ >> posi) & mask;
		const UINT32 set2 = (v_ >> posj) & mask;
		type_t xorVal = (set1 ^ set2);
		xorVal = (xorVal << posi) | (xorVal << posj);
		v_ ^= xorVal;		
	}
	
	template <typename T>
	void exchange(const T i, const T j, const T k) {
		// y[i] <--= y[k]
		// y[j] <--= y[i]
		// y[k] <--= y[j]
		const int posi = i*extent;
		const int posj = j*extent;
		const int posk = k*extent;
		const UINT32 mask = gmask0();
		const UINT32 seti = (v_ >> posi) & mask;
		const UINT32 setj = (v_ >> posj) & mask;
		const UINT32 setk = (v_ >> posk) & mask;
		type_t xorVal = (seti ^ setj);
		xorVal = (xorVal << posi) | (xorVal << posj);
		v_ ^= xorVal;
		xorVal = (setj ^ setk);
		xorVal = (xorVal << posi) | (xorVal << posk);
		v_ ^= xorVal;
	}
	
    bool operator<(const bSYT<type_t, extent>& other) const {
        return (v_<other.v_ ? true : false);
    }
    
    bool operator<=(const bSYT<type_t, extent>& other) const {
        return (v_<=other.v_ ? true : false);
    }
	
    bool operator>(const bSYT<type_t, extent>& other) const {
        return (v_>other.v_ ? true : false);
    }
    
    bool operator>=(const bSYT<type_t, extent>& other) const {
        return (v_>=other.v_ ? true : false);
    }
	
    bool operator==(const bSYT<type_t, extent>& other) const {
        return (v_==other.v_ ? true : false);
    }
	
	template <typename T>
    int operator[](const T i) const
    {
        return get(i);
    }

    void print(const unsigned int n, std::ostream& os=std::cout) const {
        os << "[" << std::flush;
        for (unsigned int i=0; i<n-1; ++i) {
            os << this->get(i) << ", " << std::flush;
        }
        os << this->get(n-1) << "]" << std::endl;
    }

protected:
	
	static UINT32 gmask0() {
		UINT32 mask = (1ULL<<extent)-1;
		return mask;
	}
	
    type_t v_;

};

#ifndef SG_BSYT_NBITS
#error "SG_BSYT_NBITS is not defined. Aborting."
#endif

#ifndef SG_BSYT_CONTAINER_SIZE
#error "SG_BSYT_CONTAINER_SIZE is not defined. Aborting."
#endif

#if SG_BSYT_CONTAINER_SIZE == 4
typedef bSYT<UINT32, SG_BSYT_NBITS> tbSYT;
#elif SG_BSYT_CONTAINER_SIZE == 8
typedef bSYT<UINT64, SG_BSYT_NBITS> tbSYT;
#endif

} // namespace sun

#endif
