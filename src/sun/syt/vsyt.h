// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_VSYT_H
#define SUN_VSYT_H

#include <iostream>
#include <vector>
#include <initializer_list>
#include <cstdint>

namespace sun {

template<class type_t>
class vSYT
{

public:
	
	vSYT() {y_.resize(0);};
    vSYT(const unsigned int n) {y_.resize(n);};
    
    template<typename T>
    vSYT(const std::vector<T> &y) {
		y_.resize(y.size());
		for (size_t i = 0; i < y.size(); ++i) {
			y_[i] = static_cast<type_t>(y[i]);
		}
	}
	
	template<typename T>
	vSYT(const std::initializer_list<T> y) {
		for (const auto & el : y) {
			y_.push_back(static_cast<type_t>(el));
		}
	}
	
	std::vector<type_t> container() const {return y_;};
	unsigned int n() const {return y_.size();};
	
	template<typename T>
    int get(const T i) const {
        return (int)y_[i];
    }
	
    template <typename T, typename U>
    void set(const T i, const U value) {
        y_[i] = static_cast<type_t>(value);
    }
	
	template <typename T>
	void exchange(const T i, const T j) {
		std::swap(y_[i], y_[j]);
	}
	
	template <typename T>
	void exchange(const T i, const T j, const T k) {
		// y[i] --> y[j]
		// y[j] --> y[k]
		// y[k] --> y[i]
		std::swap(y_[i], y_[j]);
		std::swap(y_[i], y_[k]);
	}

    bool operator<(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
        for (unsigned int i=y_.size()-1; i>0; --i) {
			if (y_[i]<other.y_[i]) {
				return true;
			} else if (y_[i]>other.y_[i]) {
				return false;
			}
		}
		return false;
    }
    
    bool operator<=(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
        for (unsigned int i=y_.size()-1; i>0; --i) {
			if (y_[i]<other.y_[i]) {
				return true;
			} else if (y_[i]>other.y_[i]) {
				return false;
			}
		}
		return true;
    }

    bool operator>(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
		for (unsigned int i=y_.size()-1; i>0; --i) {
			if (y_[i]>other.y_[i]) {
				return true;
			} else if (y_[i]<other.y_[i]) {
				return false;
			}
		}
		return false;
    }
    
    bool operator>=(const vSYT<type_t>& other) const {
		// this is the OPPOSITE of the last letter order sequence
		for (unsigned int i=y_.size()-1; i>0; --i) {
			if (y_[i]>other.y_[i]) {
				return true;
			} else if (y_[i]<other.y_[i]) {
				return false;
			}
		}
		return true;
    }

    bool operator==(const vSYT<type_t>& other) const {
        for (unsigned int i=1; i<y_.size(); ++i) {
			if (y_[i]!=other.y_[i]) {
				return false;
			}
		}
        return true;
    }
	
	template <typename T>
    int operator[](const T i) const
    {
        return (int)y_[i];
    }

    void print(const unsigned int Ns, std::ostream& os=std::cout) const {
		os << "[" << std::flush;
        for (unsigned int i=0; i<Ns-1; ++i) {
            os << y_[i] << ", " << std::flush;
        }
        os << y_[Ns-1] << "]" << std::endl;
    }

protected:

    std::vector<type_t> y_;
    
};

typedef vSYT<int8_t> Int8vSYT;
typedef vSYT<int16_t> Int16vSYT;
typedef vSYT<int32_t> Int32vSYT;

}

#endif
