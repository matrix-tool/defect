# (restricted) Defect of Unitary Matrix

## description
**Defect** of a unitary matrix `U` gives an upper bound for the dimensionality of a smooth orbit of inequivalent unitary matrices stemming from `U`. Roughly speaking, two matrices are said to be **equivalent** if one can be transformed into the other via special unitary operations. In particular it may be used as a one-way criterion
```
if defect of a matrix is zero then this matrix is an isolated point
```
which means that in the neighbourhood of a given matrix there are no other inequivalent matrices. This algebraic tool proves to be very useful when searching for potential candidates for families in the set of [complex Hadamard matrices](http://chaos.if.uj.edu.pl/~karol/hadamard/).

One can also define **restricted** defect for a unitary matrix `U` subjected to additional constraints:
- `U` is hermitian
- `U` has constant (and real) diagonal
- `U` has constant off-diagonal moduli

Motivation behind considering **restricted defect** lies in the idea of generalised quantum measurement [`POVM`](https://en.wikipedia.org/wiki/POVM) and associated Gram matrices. It turns out that Gram matrix can be expressed as a unitary and hermitian object with constant diagonal, so that the defect is applicable to it. Given a `POVM` with predefined geometrical structure - either [`SIC-POVM`](https://en.wikipedia.org/wiki/SIC-POVM) or [`MUB`](https://en.wikipedia.org/wiki/Mutually_unbiased_bases), one can ask if it is possible
to extend corresponding objects so that a smooth family of quantum measurements exists.

## implementation
Defect can be calculated in many ways. Let `R` be a matrix of a special system of linear equations associated with matrix `U` [1]. Here we present three possible implementations, appropriately named:
- 'R' - defect of `U` as the rank of `R`
- 'S' - defect of `U` as a function of non-zero singular values of `R`
- 'T' - defect of `U` as the dimension of the image of a tangent space to the manifold of unitaries under a certain tangent map... [2]

Methods 'R' and 'T' work for matrices given with the highest possible numerical precision. If matrix `U` is provided only in approximate form one can try to use the 'S' method. However, a special attention is needed when setting `SV_TOLERANCE` - a kind of threshold to distinguish between "zero" and "non-zero" singular values of the matrix `R`!

## usage
```
>> defect_u(U, [, METHOD [, SV_TOLERANCE]])

>> defect_u(U)             % implicit call of 'R' method (default)
>> defect_u(U, 'R')        % explicit call of 'R' method
>> defect_u(U, 'S')        % SVD with default SV_TOLERANCE
>> defect_u(U, 'S', 1e-12) % SVD with custom SV_TOLERANCE
>> defect_u(U, 'T')        % method of tangent spaces...
>> defect_u(U, 'R', 1e-12) % SV_TOLERANCE is ignored with 'R'
>> defect_u(U, 'T', 1e-12) % SV_TOLERANCE is ignored with 'T'
```

`defect_h.m` works similarly...

Mathematica users can find two notebooks: `defect_u.nb.txt` and `defect_h.nb`.

## references
- <sup>[1]</sup> W. Bruzda, D. Goyeneche, K. &#379;yczkowski, *"[Quantum measurements with prescribed symmetry](https://arxiv.org/abs/1704.04609 "arXiv")"*, Phys. Rev. A 96, 022105 (2017)
- <sup>[2]</sup> W. Tadej, K. &#379;yczkowski, *"[Defect of a Unitary Matrix](https://arxiv.org/abs/math/0702510 "arXiv")"*, Lin. Alg. Appl., 429, pp. 447-481 (2008)
