# Akima interpolation
# should be C1 differentiable

function akima_interpolation(sample::Real, x::Vector{Float64}, y::Vector{Float64})
    
    i = find_closest(sample, x) 
    
    m_i =  get_m(i, x, y) 
    s_i = get_s(i, x, y) # P'(x[i])
    s_i₊ = get_s(i + 1, x, y) # P'(x[i+1])

    a = y[i] # P(x[i])
    b = s_i 
    c = (3m_i - 2s_i - s_i₊) / (x[i+1] - x[i]) 
    d = (s_i + s_i₊ - 2m_i) / (x[i+1] - x[i])^2 

    sample_y = a + b * (sample - x[i]) + c * (sample - x[i])^2 + d * (sample - x[i])^3
    return sample_y

end

function get_m(i::Int, x::Vector{Float64}, y::Vector{Float64})

    m_i =  (y[i+1] - y[i]) / (x[i+1] - x[i]) 
    return m_i

end

function get_s(i::Int, x::Vector{Float64}, y::Vector{Float64})

    # more slope calculations
    if i == 1

        s = get_m(1, x, y)

    elseif i == 2

        m1 = get_m(1, x, y) # i₋ = 1
        m2 = get_m(2, x, y) # i = 2
        s = (m1 + m2) / 2

    elseif i == lastindex(x) - 1

        m1 = get_m(i-1, x, y)
        m2 = get_m(i, x, y)
        s = (m1 + m2) / 2

    elseif i == lastindex(x)

        s = get_m(i-1, x, y)

    else

        m_i₊ = get_m(i+1, x, y)
        m_i = get_m(i, x, y)
        m_i₋ = get_m(i-1, x, y)
        m_i₋₂ = get_m(i-2, x, y)

        denominator = abs(m_i₊ - m_i) + abs(m_i₋ - m_i₋₂)
        
        if denominator == 0

            s = (m_i₋ + m_i)/2

        else # weighted average

            numerator = abs(m_i₊ - m_i) * m_i₋ + abs(m_i₋ - m_i₋₂) * m_i
            s = numerator / denominator

        end

    end

    return s

end