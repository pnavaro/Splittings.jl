# Splittings.jl Documentation

Operators splitting package to solve equations of the form 

```\\frac{dU}{dt} = (T+V)U ``,

where ``T`` and  ``V`` are two differential operators by solving successively the simpler
equations

The abstract type `OperatorSplitting` implements the composition form of
different kinds of composition methods defined by their coefficients.

The application of an operator splitting method to a concrete problem is done
by extending this type containing on the one hand the data on which the operators act
and a specific implementation of the two operators

## Examples of applications are provided for

 - The linear pendulum problem.
 - The Vlasov equation with constant coefficients advection field.
 - The non linear Vlasov-Poisson equations in cartesian coordinates.

<b> References </b>

E. Hairer, C. Lubich, G. Wanner, Geometrical numerical integration, Springer 2006

This code is translated from a Fortran code written by Eric Sonnendrucker (Max-Planck-Institut für Plasmaphysik)
and Michel Mehrenberger (Aix-Marseille Université).

```@contents
```

```@meta
CurrentModule = Splittings
```

- [Semi-Lagrangian method](@ref)

## Examples

  * [Vlasov Poisson](@ref)
  * [Landau Damping](@ref)

## Functions

```@docs
compute_rho( meshv, f)
```

```@docs
compute_e( meshx, rho)
```

## Index

```@index
```
