
"""
    find_closest(x::Real, list::Vector{T}; round::Symbol=:floor)
    (find_closest(x_samples::Vector{T}, list::Vector{T}; round::Symbol=:floor))

primary method get the index of the item in the list that is closest to x. The list must be 
monotonically increasing and x must be within the bounds of the list.
(secondary method takes a list of x_samples and runs find_closest for each of them, returning
a list of indices i.)

# Arguments
- x::Real : the number you're trying to find 
- list::Vector{T} : the monotonically increasing list you're searching through
- round::Symbol=:floor : gives the option of which index to round to in the event that x falls between two list values
                         - :floor (closest index below)
                         - :ceil, :ceiling (closest index above)
"""

function find_closest(x::Real, list::Vector{T}; round::Symbol=:floor) where T <: Real

    @assert round in [:floor, :ceil, :ceiling]

    closest_i = 1

    for i in 1:lastindex(list)

        if list[i] > x

            closest_i = i-1
            return maximum([1, closest_i])

        end

    end

    return lastindex(list) - 1

end

function find_closest(x_samples::Vector{T}, list::Vector{T}; round::Symbol=:floor) where T <: Real

    indices = [find_closest(x_samples[i],list,round=round) for i in 1:length(x_samples)]
    return indices

end