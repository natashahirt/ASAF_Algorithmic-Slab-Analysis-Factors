data = XLSX.readdata(joinpath(@__DIR__, "data/aisc-shapes-database-v15.0.xlsx"), "Database v15.0!A2:ED2041")

const Wrange = 1:283
const Crange = 352:383
const Lrange = 424:560
const WTrange = 561:843
const LLrange = 886:1524
const HSSRectrange = 1525:1912
const HSSRoundrange = 1913:2040

# AISC v15.0 Column Mapping
# CH is column 86. Our range starts at A (1).
const OFFSET = 85

# column index dictionary
colDict = Dict(
    :name => 1 + OFFSET,          # EDI_Std_Nominal (Metric EDI Name)
    :name_imperial => 2,          # AISC_Manual_Label (Imperial Name) - Column B
    :W => 2 + OFFSET,
    :A => 3 + OFFSET,
    :d => 4 + OFFSET,
    :b => 12 + OFFSET,
    :bf => 9 + OFFSET,
    :tw => 14 + OFFSET,
    :tf => 17 + OFFSET,
    :t => 19 + OFFSET,
    :Ix => 36 + OFFSET,
    :Zx => 37 + OFFSET,
    :Sx => 38 + OFFSET,
    :rx => 39 + OFFSET,
    :Iy => 40 + OFFSET,
    :Zy => 41 + OFFSET,
    :Sy => 42 + OFFSET,
    :ry => 43 + OFFSET,
    :J => 47 + OFFSET,
    :Cw => 48 + OFFSET,
    :Ht => 6 + OFFSET,
    :B => 11 + OFFSET,
    :tHSS => 21 + OFFSET,
    :OD => 8 + OFFSET
    )

names = data[:, colDict[:name]]

const field_units = Dict(
    :name => nothing,
    :name_imperial => nothing,
    :A => u"mm^2",
    :d => u"mm",
    :b => u"mm",
    :bf => u"mm",
    :tw => u"mm",
    :tf => u"mm",
    :t => u"mm",
    :Ix => 1e6 * u"mm^4",
    :Zx => 1e3 * u"mm^3",
    :Sx => 1e3 * u"mm^3",
    :rx => u"mm",
    :Iy => 1e6 * u"mm^4",
    :Zy => 1e3 * u"mm^3",
    :Sy => 1e3 * u"mm^3",
    :ry => u"mm",
    :J => 1e3 * u"mm^4",
    :Cw => 1e9 * u"mm^6",
    :Ht => u"mm",
    :B => u"mm",
    :tHSS => u"mm",
    :OD => u"mm"
)

const Wfields = [:name, :name_imperial, :A, :d, :bf, :tw, :tf, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry, :J, :Cw]
const Cfields = [:name, :name_imperial, :A, :d, :bf, :tw, :tf, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry, :J, :Cw]
const Lfields = [:name, :name_imperial, :A, :b, :d, :t, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry, :J, :Cw]
const LLfields = [:name, :name_imperial, :A, :b, :d, :t, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry]
const HSSRectfields = [:name, :name_imperial, :A, :Ht, :B, :tHSS, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry, :J]
const HSSRoundfields = [:name, :name_imperial, :A, :OD, :tHSS, :Ix, :Zx, :Sx, :rx, :Iy, :Zy, :Sy, :ry, :J]

macro SteelSections(typeName, superType, fields, range)
    return quote
        struct $typeName <: $superType
            name::String
            name_imperial::String
            # Generate fields from the field list
            $(Expr(:block, [Expr(:(::), f, :Any) for f in eval(fields)[3:end]]...))

            function $(typeName)(name::String)
                irow = findfirst(names .== name)
                if isnothing(irow)
                    imperial_names = data[:, colDict[:name_imperial]]
                    irow = findfirst(imperial_names .== name)
                end
                @assert !isnothing(irow) "Section name '$name' not found in Metric or Imperial records."
                
                vals = Any[data[irow, colDict[f]] for f in $fields]
                
                # Apply Unitful units and multipliers
                for i in 3:length(vals)
                    f = $(fields)[i]
                    u = get(field_units, f, nothing)
                    if !isnothing(u)
                        vals[i] *= u
                    end
                end
                
                return new(vals...)
            end
        end
    end |> esc
end
