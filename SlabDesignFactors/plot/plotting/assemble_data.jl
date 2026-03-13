"""
    assemble_data(files::Vector{String}; category_names::Vector{String}=[]) -> DataFrame

Assembles data from multiple CSV files into a single DataFrame with calculated metrics and standardized formatting.

# Arguments
- `files::Vector{String}`: Vector of paths to CSV files containing slab design data
- `category_names::Vector{String}=[]`: Optional vector of category names for each file. If empty, uses filenames as categories.

# Returns
- `DataFrame`: Combined DataFrame with calculated metrics including:
  - Embodied carbon metrics (steel_ec, concrete_ec, rebar_ec, total_ec)
  - Geometric properties (symbol, rotation, vector coordinates)
  - Row/column identifiers
  - Original data from CSV files

The function processes each file, standardizes the data format, calculates derived metrics,
and applies consistent geometric properties based on slab type. It filters out rows with zero area
and sorts by row/column position.
"""
function assemble_data(files::Vector{String}; category_names::Vector{String}=String[])

    if isempty(category_names)
        category_names = [replace(basename(file), ".csv" => "") for file in files]
    else
        @assert length(category_names) == length(files) "You need the same number of category names as files"
    end

    df = DataFrame(category=String[], name=String[], area=Float64[], steel_norm=Float64[], concrete_norm=Float64[], rebar_norm=Float64[], max_depth=Float64[], slab_type=Symbol[], slab_sizer=Symbol[], beam_sizer=Symbol[], collinear=Bool[], symbol=Symbol[], rotation=Float64[], vector_1d_x=Float64[], vector_1d_y=Float64[], sections=Any[], ids=Any[], unique_sections=Int64[])

    for (i,filename) in enumerate(files)

        df_slab = DataFrame(CSV.File(filename))
        df_slab.symbol .= :circle
        df_slab.rotation .= 0.
        df_slab.category .= category_names[i]
        if !hasproperty(df_slab, :unique_sections)
            df_slab.unique_sections .= 0
        end

        df = vcat(df, df_slab)
    
    end

    df.steel_ec = df.steel_norm .* ECC_STEEL
    df.concrete_ec = df.concrete_norm .* ECC_CONCRETE
    df.rebar_ec = df.rebar_norm .* ECC_REBAR
    df.slab_ec = df.concrete_ec + df.rebar_ec
    df.total_ec = df.steel_ec + df.concrete_ec + df.rebar_ec
    df.row .= 0
    df.col .= 0
    df.rowcol .= ""

    for i in 1:lastindex(df.name)

        row = df[i,:]

        if row.slab_type == "isotropic"

            row.symbol = :star8
            row.rotation = 0.
            row.vector_1d_x = 0.
            row.vector_1d_y = 0.

        elseif row.slab_type == "orth_biaxial"

            row.symbol = :cross
            row.rotation = get_vector_1d_angle([row.vector_1d_x,row.vector_1d_y])

        elseif row.slab_type == "uniaxial"

            row.symbol = :hline
            row.rotation = get_vector_1d_angle([row.vector_1d_x,row.vector_1d_y])

        end

        # Check if name has format "rXcY" where X and Y are numbers
        if occursin(r"^[a-zA-Z]?\d+[a-zA-Z]?\d+$", row.name)
            split_name = split(row.name, r"(?<=\d)(?=\D)|(?<=\D)(?=\d)") # \d is decimal digit, \D is nondigit characters
            row.row = parse(Int,split_name[2])
            row.col = parse(Int,split_name[4])
        else
            row.row = 0
            row.col = 0
        end 

        # bump one up for the grid
        if contains(row.name, "x" )&& contains(row.name, "y")
            row.row += 1
            row.col += 1
        end

        row.rowcol = "$(row.category)[$(row.row),$(row.col)]"
        row.ids = split(row.ids, ",")

    end

    println("values: ",length(df.name))

    df = filter(row -> row.area != 0 && !isnan(row.area), df)
    sort!(df, [:row, :col])

    GC.gc()

    return df

end

 function assemble_data(file::String)
    if isdir(file) || !endswith(file, ".csv")
        # Get all CSV files in directory
        files = filter(f -> endswith(f, ".csv"), readdir(file))
        
        # Call csv_2_df with the full file paths
        files = String[joinpath(file, f) for f in files if !contains(f, "max_depth") && !contains(f, "fix_params")]
        return assemble_data(files)
    else
        # Single file case
        return assemble_data(String[(file)])
    end
end

function assemble_data(dfs::Vector{DataFrame}, category_name::String, boolean::Vector{Bool})

    for (i,df) in enumerate(dfs)
        df[!, category_name] .= boolean[i]
    end

    full_df = vcat(dfs...)

    return full_df

end