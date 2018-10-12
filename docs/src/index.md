# Splittings.jl Documentation

Operators splitting package to solve equations of the form 

`` dU/dt = (T+V)U ``,

where ``T`` and  ``V`` are two differential operators by solving successively the simpler
equations

The application of an operator splitting method to a concrete problem is done
by using Julia macros.

## Examples of applications are provided for

 - The linear pendulum problem.
 - The Vlasov equation with constant coefficients advection field.
 - The non linear Vlasov-Poisson equations in cartesian coordinates.

*References*

E. Hairer, C. Lubich, G. Wanner, Geometrical numerical integration, Springer 2006

This code is derived from Fortran and Python codes written by 

    - Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences).
    - Michel Mehrenberger  (Aix-Marseille Université).
    - Eric Sonnendrucker   (Max-Planck-Institut für Plasmaphysik).

```@contents
```

```@meta
CurrentModule = Splittings
```

- [Semi-Lagrangian method](@ref)

## Examples

  * [Vlasov-Poisson](@ref)
  * [Vlasov-Ampere](@ref)
  * [Bump On Tail](@ref)

## Functions

```@autodocs
Modules = [Splittings]
Order   = [:advection!, :UniformMesh]
```

```@docs
compute_rho( meshv, f)
```

```@docs
compute_e( meshx, rho)
```

## Index

```@index
```
