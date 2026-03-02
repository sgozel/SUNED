// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_VSYT_H
#define SUN_VSYT_H

#include <iostream>
#include <vector>
#include <cstdint>

namespace sun {

template<class type_t>
class vSYT
{

public:
	
	vSYT() {y = std::vector<type_t>(0);};
    vSYT(const unsigned int n) {y = std::vector<type_t>(n);};
	
	template<typename T>
    int get(const T i) const {
        return (int)y[i];
    }
	
    template <typename T, typename U>
    void set(const T i, const U value) {
        y[i] = static_cast<type_t>(value);
    }
	
	template <typename T>
	void exchange(const T i, const T j) {
		std::swap(y[i], y[j]);
	}
	
	template <typename T>
	void exchange(const T i, const T j, const T k) {
		// y[i] --> y[j]
		// y[j] --> y[k]
		// y[k] --> y[i]
		std::swap(y[i], y[j]);
		std::swap(y[i], y[k]);
	}

    bool operator<(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
        for (unsigned int i=y.size()-1; i>0; --i) {
			if (y[i]<other.y[i]) {
				return true;
			} else if (y[i]>other.y[i]) {
				return false;
			}
		}
		return false;
    }

    bool operator>(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
		for (unsigned int i=y.size()-1; i>0; --i) {
			if (y[i]>other.y[i]) {
				return true;
			} else if (y[i]<other.y[i]) {
				return false;
			}
		}
		return false;
    }

    bool operator==(const vSYT<type_t>& other) const {
        for (unsigned int i=1; i<y.size(); ++i) {
			if (y[i]!=other.y[i]) {
				return false;
			}
		}
        return true;
    }
	
	template <typename T>
    int operator[](const T i) const
    {
        return (int)y[i];
    }

    void print(const unsigned int Ns, std::ostream& os=std::cout) const {
		os << "[" << std::flush;
        for (unsigned int i=0; i<Ns-1; ++i) {
            os << y[i] << ", " << std::flush;
        }
        os << y[Ns-1] << "]" << std::endl;
    }

protected:

    std::vector<type_t> y;
    
};

typedef vSYT<int8_t> Int8vSYT;
typedef vSYT<int16_t> Int16vSYT;
typedef vSYT<int32_t> Int32vSYT;

}

#endif
