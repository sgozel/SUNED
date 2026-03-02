// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_PERMUTATION_H
#define SUN_PERMUTATION_H

#include <iostream>
#include <initializer_list>
#include <vector>


namespace sun {

class Transposition;
class AdjacentTransposition;

class Cycle
{
public:
    Cycle();
    Cycle(const std::initializer_list<int> &);
    Cycle(const std::vector<int> &);

    std::vector<Transposition> toTransposition() const;
    virtual std::vector<AdjacentTransposition> toAdjacentTransposition() const;
	
    friend std::ostream& operator<<(std::ostream&, const Cycle&);

public:
    std::vector<int> cycle_;
};


class Transposition : public Cycle
{
public:
    Transposition();
    Transposition(const int i, const int j);
    Transposition(const std::initializer_list<int> &);
    Transposition(const std::vector<int> &);
    
    int geti() const { return cycle_[0]; };
    int getj() const { return cycle_[1]; };

    virtual std::vector<AdjacentTransposition> toAdjacentTransposition() const;
};


class AdjacentTransposition : public Transposition
{
public:
    AdjacentTransposition();
    AdjacentTransposition(const int k); // Transposition (k, k+1)

    int getk() const {return cycle_[0];};

    virtual std::vector<AdjacentTransposition> toAdjacentTransposition() const;
};


template<typename T, typename = std::enable_if_t<std::is_base_of_v<Cycle, T>>>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& cycles)
{
    for (const Cycle& cycle : cycles) {
        os << cycle;
    }
    return os;
}

using Permutation = std::vector<Cycle>;

std::vector<AdjacentTransposition> toAdjacentTransposition(const Permutation&);

void print(const Permutation & p);

void print(const std::vector<AdjacentTransposition> & adjaTranspos);

} // namespace sun

#endif
