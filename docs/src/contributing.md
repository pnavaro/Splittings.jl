# How to Contribute

Here's an outline of the workflow you should use if you want to make 
contributions to Splittings.

1. Fork Splittings
2. Make a new branch on your fork, named after whatever changes you'll be making
3. Apply your code changes to the branch on your fork
4. When you're done, submit a PR to Splittings to merge your fork into 
   Splittings's master branch.


# Add a new feature 

Adding new functions to a Julia package is not easy because you must reload
the package to test your code. The best way is to run julia in
`Splittings` directory in pkg mode and activate the package:

```julia
julia> pwd()
"/Users/navaro/.julia/dev/Splittings"

(v1.0) pkg> activate .

(Splittings) pkg>
```

Write your program in `test` directory and add an `include` line in the
file `runtests.jl`. Just type test in the Julia console.

```julia
(Splittings) pkg> test
   Testing Splittings
 Resolving package versions...
```

When it is finished you can move your functions in `src` directory and just keep
the test.

