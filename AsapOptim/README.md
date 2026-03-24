![](figures/gradients-axo.png)

# AsapOptim

High performance *general* structural optimization in the [Asap.jl](https://github.com/keithjlee/Asap) environment. It provides the following:

- A set of data structures for controllable, extendable, and composable definitions of design variables and objective functions
- `CoupledVariable`s that reduce the dimensionality of design problems
- A complete set of custom chain rules for extremely high-performing automatic differentiation during gradient calculations
- A set of differentiable definitions of common structural objectives

This is highly experimental and prone to breaking changes. Currently unstable for Julia versions >= 1.9. 