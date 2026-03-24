abstract type LoadVariableProperty end

"""
    LoadVariable <: IndependentVariable

A variable representing a vertical load applied to an element.

```
LoadVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :LineLoad)
LoadVariable(element::Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :LineLoad)
LoadVariable(element::Element, lowerbound::Float64, upperbound::Float64, property::Symbol = :LineLoad)
```
"""

mutable struct LoadVariable{T<:LoadVariableProperty} <: IndependentVariable
    i::Int64
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Float64
end

const property_to_load_type = Dict(
    :LineLoad => LineLoad
)

function LoadVariable(element::Asap.Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :LineLoad)

    T = property_to_load_type[property]
    
    return LoadVariable{T}(element.elementID, value, lowerbound, upperbound, 0)

end