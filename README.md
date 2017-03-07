# (restricted) Defect of Unitary Matrix

## description
**Defect** of a unitary matrix `U` gives an upper bound for the dimensionality of a smooth orbit of inequivalent unitary matrices stemming from `U`. In particular one can use it as a one-way criterion
```
if defect of a matrix is zero then this matrix is an isolated point
```
which means that in the neighbourhood of a given matrix there are no other inequivalent matrices. Roughly speaking, two matrices are said to be **equivalent** is one can be transformed into the other via unitary operations ("rotations"). This algebraic tool proves to be very useful when searching for potential candidates for families in the set of [complex Hadamard matrices](http://chaos.if.uj.edu.pl).

One can also define **restricted** defect for a unitary matrix `U` subjected to additional constraints:
- `U` is hermitian
- `U` has constant (and real) diagonal
- `U` has constant off-diagonal moduli

Motivation behind considering **restricted defect** is the idea of generalised quantum measurement and associated Gram matrices. It turns out that Gram matrix can be expressed as a unitary and hermitian object with constant diagonal, so that the defect is applicable to it. Given a POVM with predefined geometrical structure - either `SIC-POVM` or `MUB`, one can ask if it is possible
to extend corresponding objects so that a smooth family of quantum measurements exists. In many cases the answer is affirmative.

## implementation
Defect can be calculated in many ways. Let `R` be a matrix of a special system of linear equations associated with matrix `U` [1]. Here we present three possible implementations, appropriately named:
- 'R' - defect of `U` as the rank of `R`
- 'S' - defect of `U` as a function of non-zero singular values of `R`
- 'T' - defect of `U` as the dimension of the image of a tangent space to the manifold of unitaries under a certain tangent map... [2]

Methods 'R' and 'T' work for matrices given with the highest possible numerical precision. If matrix `U` is provided only in approximate form one can try to use the 'S' method. However, a special attention is needed when setting `SV_TOLERANCE` - a kind of threshold to distinguish between "zero" and "non-zero" singular values of the matrix `R`!

## usage
```
>> defect(U, [, METHOD [, SV_TOLERANCE]])

>> defect(U)             % implicit call of 'R' method (default)
>> defect(U, 'R')        % explicit call of 'R' method
>> defect(U, 'S')        % SVD with default SV_TOLERANCE
>> defect(U, 'S', 1e-12) % SVD with custom SV_TOLERANCE
>> defect(U, 'T')        % method of tangent spaces...
>> defect(U, 'R', 1e-12) % SV_TOLERANCE is ignored with 'R'
>> defect(U, 'T', 1e-12) % SV_TOLERANCE is ignored with 'T'
```

The main script `defect.m` comprises of two auxiliary scripts which can be used independently. They are:
- `defect_u.m` - defect of a unitary matrix
- `defect_h.m` - defect of a unitary matrix subjected to additional constraints

## references
- <sup>[1]</sup> W. Bruzda, D. Goyeneche, K. Życzkowski, *"Quantum measurements with prescribed symmetry"*, preprint (2017)
- <sup>[2]</sup> W. Tadej, K. Życzkowski, *"[Defect of a Unitary Matrix](https://arxiv.org/abs/math/0702510 "arXiv")"*, Linear Algebra and its Applications, 429, pp. 447-481 (2008)
