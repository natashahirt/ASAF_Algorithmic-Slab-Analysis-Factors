"""
    to_imperial(self::InternalForces) -> InternalForces

Converts the internal forces from metric units to imperial units.

# Arguments
- `self::InternalForces`: An instance of `InternalForces` containing metric units.

# Returns
- `InternalForces`: A new instance with forces converted to imperial units.
"""
function to_imperial(self::InternalForces)
    # Extract properties from the input object
    element = self.element
    resolution = self.resolution

    # Convert units
    x = self.x .* (3.280840 * 12) # meters to inches
    P = self.P .* 0.224809 # kN to kips
    My = self.My .* (0.737562 * 12) # kNm to kip-in
    Mz = self.Mz .* (0.737562 * 12) # kNm to kip-in
    Vy = self.Vy .* 0.224809 # kN to kips
    Vz = self.Vz .* 0.224809 # kN to kips    

    # Create a new InternalForces object with converted units
    imperial_forces = InternalForces(element, resolution, x, P, My, Vy, Mz, Vz)

    return imperial_forces
end

"""
    to_metric(self::InternalForces) -> InternalForces

Converts the internal forces from imperial units to metric units.

# Arguments
- `self::InternalForces`: An instance of `InternalForces` containing imperial units.

# Returns
- `InternalForces`: A new instance with forces converted to metric units.
"""
function to_metric(self::InternalForces)
    # Extract properties from the input object
    element = self.element
    resolution = self.resolution

    # Convert units
    x = self.x ./ (3.280840 * 12) # inches to meters
    P = self.P ./ 0.224809 # kips to kN
    My = self.My ./ (0.737562 * 12) # kip-in to kNm
    Mz = self.Mz ./ (0.737562 * 12) # kip-in to kNm
    Vy = self.Vy ./ 0.224809 # kips to kN
    Vz = self.Vz ./ 0.224809 # kips to kN

    # Create a new InternalForces object with converted units
    metric_forces = InternalForces(element, resolution, x, P, My, Vy, Mz, Vz)

    return metric_forces
end
