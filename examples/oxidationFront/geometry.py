import viennals as vls

def get_geometry_domains(params):
    """
    Generates the LevelSets and MaterialMap using the provided parameters.
    Expects params to contain:
      - gridDelta, xExtent, substrateHeight, maskHeight, holeRadius
    
    Returns: (list_of_levelsets, material_map, bounds, boundary_conditions, grid_delta)
    """
    print("--- Generating Oxidation Geometry (ViennaLS) ---")

    # Unpack parameters
    grid_delta = params["gridDelta"]
    x_extent = params["xExtent"]
    substrate_height = params["substrateHeight"]
    mask_height = params["maskHeight"]
    hole_radius = params["holeRadius"]

    # Define Bounds
    # Domain height covers substrate + mask + buffer
    y_max = substrate_height + mask_height + grid_delta
    bounds = [-x_extent/2.0, x_extent/2.0, 0.0, y_max]
    
    bc = [
        vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        vls.BoundaryConditionEnum.INFINITE_BOUNDARY
    ]

    # --- Layer 0: Bottom Plane (y=0) ---
    bottom = vls.Domain(bounds, bc, grid_delta)
    vls.MakeGeometry(bottom, vls.Plane(origin=[0, 0], normal=[0, 1])).apply()

    # --- Layer 1: Substrate Surface (y=substrateHeight) ---
    substrate = vls.Domain(bounds, bc, grid_delta)
    vls.MakeGeometry(substrate, vls.Plane(origin=[0, substrate_height], normal=[0, 1])).apply()
    
    # Union Bottom and Substrate
    vls.BooleanOperation(substrate, bottom, vls.BooleanOperationEnum.UNION).apply()

    # --- Layer 2: Mask (y=substrateHeight + maskHeight) ---
    mask = vls.Domain(bounds, bc, grid_delta)
    vls.MakeGeometry(mask, vls.Plane(origin=[0, substrate_height + mask_height], 
                                     normal=[0, 1])).apply()
    
    # Hole
    hole = vls.Domain(bounds, bc, grid_delta)
    min_p = [-hole_radius, substrate_height - grid_delta]
    max_p = [hole_radius, substrate_height + mask_height + grid_delta]
    vls.MakeGeometry(hole, vls.Box(minPoint=min_p, maxPoint=max_p)).apply()
    
    # Cut Hole
    vls.BooleanOperation(mask, hole, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()
    
    # Union Mask with Substrate
    vls.BooleanOperation(mask, substrate, vls.BooleanOperationEnum.UNION).apply()

    # --- Prepare Output ---
    ls_list = [bottom, substrate, mask]
    
    mat_map = vls.MaterialMap()
    mat_map.insertNextMaterial(0) # Substrate
    mat_map.insertNextMaterial(0) # Substrate (redundant layer)
    mat_map.insertNextMaterial(1) # Mask

    print("Geometry generated.")
    return ls_list, mat_map, bounds, bc, grid_delta

if __name__ == "__main__":
    # Test with default values if run directly
    test_params = {
        "gridDelta": 0.85,
        "xExtent": 85.0,
        "substrateHeight": 50.0,
        "maskHeight": 10.0,
        "holeRadius": 20.0
    }
    lss, _, _, _, _ = get_geometry_domains(test_params)
    mesh = vls.Mesh()
    vls.ToSurfaceMesh(lss[1], mesh).apply()
    vls.VTKWriter(mesh, "check_substrate.vtp").apply()
    print("Written check_substrate.vtp using test params")