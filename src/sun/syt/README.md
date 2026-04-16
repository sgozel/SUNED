# SYT classes

`SUNED` has been developed to allow different storing strategies for Standard Young Tableaux (SYTs), the building block of ED with SU($N$) symmetry, provided they defined a minimum set of APIs. We present here two strategies. The first one, `vSYT`, is the most natural and straighforward approach. `bSYT` is our new cutting-edge strategy, which is both faster and lighter than `vSYT`.

## vSYT

`vsyt.h` defines a standard `v`ector strategy for storing Standard Young Tableaux (SYTs). A SYT for an irrep with $n$ boxes is stored in a `std::vector` of length $n$. The class is templated on the box-container type, but `int8_t` should be the preferred choice. The class implements all necessary operations on SYTs: get, set and 2- and 3-way exchange of boxes. It also overloads the comparison operators, following the *decreasing* order of the last letter order sequence.

## bSYT

`bsyt.h` defines our new light and efficient bit-level storing strategy for Standard Young Tableaux (SYTs) of SU($N$). It implements a SYT in a *single* word (container), whose number of bytes must be provided at compile time with the `SG_BSYT_CONTAINER_SIZE` macro. When used in the `SUNED` code, this macro is controlled through the CMake option `SUNED_BSYT_CONTAINER_SIZE`. The box positions are packed in a sequential way, ensuring extremely fast bit-level get, set and exchanges. Each box requires $`\eta := \lceil \log_2(N_r) \rceil`$ bits, where $`N_r \leq N`$ is the number of rows in the Young diagram. In practice, for all relevant systems to be numerically investigated, one actually has $N_r = N$. The `bSYT` class is templated on the container size and on the number $`\eta`$ of bits per box. The latter needs also to be provided at compile time through the `SG_BSYT_NBITS` macro, which is itself controlled by the CMake option `SUNED_BSYT_NBITS`. Having both these macro defined at compile time ensures maximal performance for operations on SYTs. A great benefit of `bSYT` compared to `vSYT` is that it implies a *natural* ordering of SYTs through their container value. In fact, it is easy to convince oneself that this order is the *decreasing* order of the last letter order sequence. As a consequence, searching across an ordered list of SYTs becomes as trivial as performant: a binary search with *single-container comparison* at each step.

In summary, the benefits of using the new `bSYT` storing strategy versus the standard `vSYT` are:
- Memory efficiency: `vSYT` uses $n$ bytes per SYT, where $n$ is the system size. `bSYT` uses a fixed number of bytes per SYT ($4$, $8$ or $16$). How this translates in practice ? Take the $n=28$ SU(4) system on the singlet irrep. Storing all 13.67 billion SYTs using `vSYT` amounts to 383 GB of memory. With `bSYT`, the entire basis is stored with 109 GB, that's a 3.5x memory gain.

- Speed: the details about the container for a `bSYT` are selected at compile time. All operations on `bSYT`s are implemented at the bit level. How this translates in practice ? Generating all 140 millions SYTs for the SU(4) singlet with 24 sites takes 11.7 seconds with `vSYT`. It takes 3.7 seconds with `bSYT`, that's a 3.2x speedup.

- Speed again: *indexing* a SYT, namely finding the index of this SYT in the ordered basis of all SYTs is one of the most important operations in light of Exact Diagonalization. For `vSYT`, this is particularly expensive, as a single comparison operator on two `vSYT`s involves a loop with up to $n-1$ "elementary" comparisons on the box positions of the two `vSYT`s. For `bSYT`, this is replaced by a single comparison on a single word. How this translates in practice ? When indexing an entire basis on a single thread, the speedup is around 6.2x with `bSYT`. The lighter memory structure of `bSYT` means it is also more suited for parallelization, before being limited by memory bandwidth.

## License

The code is licensed under GNU GPL-v3.0 as given in the file LICENSE.

## Author

Samuel Gozel
