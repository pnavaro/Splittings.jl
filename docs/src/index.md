# Splittings.jl Documentation

Operators splitting package to solve equations of the form 
```math
 \frac{dU}{dt} = (T+V)U,
```
where T and V are differential operators.

The solution on one time step can be written 
```math
U(Δt) = \mathcal{S}_{T+V} U(0). 
```
The composition algorithm consists in successive solutions of the split equations 
```math
\frac{dU}{dt} = T U 
```
and 
```math
\frac{dU}{dt} = V U 
``` 

Alternating the two reduced solution operators ``\mathcal{S}_{T}``
and ``\mathcal{S}_{V}`` with adequately chosen time increments yields arbitrary 
order in time for the full solution.

The application of an operator splitting method to a concrete problem is done
by using Julia macros:

```@docs
@Lie
@Strang
@TripleJump
@Order6
```

Examples of applications are provided for~:

 - The linear pendulum problem.
 - The Vlasov equation with constant coefficients advection field.
 - The non linear Vlasov-Poisson equations in cartesian coordinates.

This code is derived from Fortran and Python codes written by~:

 - Edwin Chacon Golcher (Institute of Physics of the Czech Academy of Sciences).
 - Michel Mehrenberger  (Aix-Marseille Université).
 - Eric Sonnendrucker   (Max-Planck-Institut für Plasmaphysik).
