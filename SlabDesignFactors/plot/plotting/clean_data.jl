# Fix the ids != sections bug when collinear was saved
# Read CSV into DataFrame first to allow modification
function clean_csv_data(file::String)

    if file[end] == '/'
        for category in ["topology", "grid", "nova"]
            filename = file * category * ".csv"
            clean_csv_data(filename)
        end
        return
    end

    df = CSV.read(file, DataFrame)

    for i in 1:nrow(df)
        sections = parse_sections(String(df[i, :sections])) 
        ids = parse_sections(String(df[i, :ids]))
        if sections != ids
            df[i, :sections] = "Any[" * join(map(x -> "\"$x\"", ids), ", ") * "]"
        end
    end

    # Write back to CSV if needed
    CSV.write(file, df)

    for row in CSV.Rows(file)
        sections = parse_sections(String(row.sections)) # Convert PosLenString to String before parsing
        ids = parse_sections(String(row.ids))
        if sections != ids
            println(row.name)
        end
    end

end