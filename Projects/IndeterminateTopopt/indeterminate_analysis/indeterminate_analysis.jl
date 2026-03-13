function calculate_slab_loads_indeterminate(self::SlabAnalysisParams)

    self.model.loads = Asap.AbstractLoad[]

    if isempty(self.raster_df)
        self = get_raster_df(self; resolution=50)
    end

    self = get_load_probabilities(self)
    self = get_raster_loads(self)

    self.load_dictionary = get_load_dictionary_by_id(self)

    println("Number of model loads: $(length(self.model.loads))")
    println("Slab area: $(round(sum(self.areas), digits=2)) m^2")
    
    return self

end
