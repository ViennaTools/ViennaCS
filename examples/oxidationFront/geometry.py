import sys
import os

import viennals as vls

def make_structure(params, mat_map, substrate_material=0, mask_material=2, ambient_material=3):
    """
    Generates the LevelSets using the provided parameters and material map.
    Matches the C++ geometry::makeStructure template function (geometry.hpp lines 55-143).

    For 2D: Creates plane-based geometry with bottom, substrate, mask with hole, and ambient.
    For 3D: Creates cylindrical geometry with substrate, mask with hole, and ambient.

    Args:
        params: Dictionary containing geometry parameters
        mat_map: vls.MaterialMap object to populate
        vls: viennals2d or viennals3d module
        substrate_material: Material ID for substrate (default: 0)
        mask_material: Material ID for mask (default: 2)
        ambient_material: Material ID for ambient (default: 3)

    Returns: list_of_levelsets
    """
    # Unpack parameters
    dimension = int(params.get("dimensions", params.get("dimension", 3)))
    vls.setDimension(dimension)
    grid_delta = params["gridDelta"]
    x_extent = params["xExtent"]
    y_extent = params.get("yExtent", x_extent)  # For 3D
    substrate_height = params["substrateHeight"]
    mask_height = params["maskHeight"]
    hole_radius = params["holeRadius"]
    ambient_height = params.get("ambientHeight", params.get("coverHeight", 0.0))

    # Define Bounds
    # Domain height covers substrate + mask + ambient + buffer
    z_max = substrate_height + mask_height + ambient_height + grid_delta

    if dimension == 2:
        bounds = [-x_extent/2.0, x_extent/2.0, 0.0, z_max]
        bc = [
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.INFINITE_BOUNDARY
        ]
    else:  # dimension == 3
        bounds = [-x_extent/2.0, x_extent/2.0, -y_extent/2.0, y_extent/2.0, 0.0, z_max]
        bc = [
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.INFINITE_BOUNDARY
        ]

    level_sets = []
    origin = [0.0] * dimension
    normal = [0.0] * dimension
    normal[-1] = 1.0  # Normal points in Z direction

    # --- Substrate ---
    if dimension == 2:
        # Layer 0: Bottom Plane (y=0)
        bottom = vls.Domain(bounds, bc, grid_delta)
        vls.MakeGeometry(bottom, vls.Plane(origin=origin, normal=normal)).apply()
        level_sets.append(bottom)
        mat_map.insertNextMaterial(substrate_material)

        # Layer 1: Substrate Surface (y=substrateHeight)
        origin[-1] = substrate_height
        substrate = vls.Domain(bounds, bc, grid_delta)
        vls.MakeGeometry(substrate, vls.Plane(origin=origin, normal=normal)).apply()
        vls.BooleanOperation(substrate, bottom, vls.BooleanOperationEnum.UNION).apply()
        level_sets.append(substrate)
        mat_map.insertNextMaterial(substrate_material)
    else:  # dimension == 3
        origin[-1] = 0.0
        substrate = vls.Domain(bounds, bc, grid_delta)
        radius = x_extent / 2.0
        vls.MakeGeometry(substrate, vls.Cylinder(origin=origin, axisDirection=normal,
                                                   height=substrate_height, radius=radius)).apply()
        level_sets.append(substrate)
        mat_map.insertNextMaterial(substrate_material)

    # --- Mask (if maskHeight > 0) ---
    if mask_height > 0:
        mask = vls.Domain(bounds, bc, grid_delta)

        if dimension == 2:
            origin[-1] = substrate_height + mask_height
            vls.MakeGeometry(mask, vls.Plane(origin=origin, normal=normal)).apply()

            # Create hole in mask using Box
            hole = vls.Domain(bounds, bc, grid_delta)
            min_p = [-hole_radius, substrate_height - grid_delta]
            max_p = [hole_radius, substrate_height + mask_height + grid_delta]
            vls.MakeGeometry(hole, vls.Box(minPoint=min_p, maxPoint=max_p)).apply()
        else:  # dimension == 3
            origin[-1] = substrate_height
            radius = x_extent / 2.0
            vls.MakeGeometry(mask, vls.Cylinder(origin=origin, axisDirection=normal,
                                                  height=mask_height, radius=radius)).apply()

            # Create hole in mask using smaller Cylinder
            hole = vls.Domain(bounds, bc, grid_delta)
            vls.MakeGeometry(hole, vls.Cylinder(origin=origin, axisDirection=normal,
                                                  height=mask_height + 2*grid_delta,
                                                  radius=hole_radius)).apply()

        # Cut hole from mask
        vls.BooleanOperation(mask, hole, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()

        # Union mask with substrate
        vls.BooleanOperation(mask, level_sets[-1], vls.BooleanOperationEnum.UNION).apply()
        level_sets.append(mask)
        mat_map.insertNextMaterial(mask_material)

    # --- Ambient (Gas phase) ---
    if ambient_height > 0:
        ambient = vls.Domain(bounds, bc, grid_delta)

        if dimension == 2:
            origin[-1] = substrate_height + mask_height + ambient_height
            vls.MakeGeometry(ambient, vls.Plane(origin=origin, normal=normal)).apply()
        else:  # dimension == 3
            origin[-1] = substrate_height
            radius = x_extent / 2.0
            vls.MakeGeometry(ambient, vls.Cylinder(origin=origin, axisDirection=normal,
                                                     height=ambient_height, radius=radius)).apply()

        # Union ambient with previous layer (following C++ addLevelSet pattern)
        if len(level_sets) > 0:
            vls.BooleanOperation(ambient, level_sets[-1], vls.BooleanOperationEnum.UNION).apply()
        level_sets.append(ambient)
        mat_map.insertNextMaterial(ambient_material)

    return level_sets

def get_geometry_domains(params, substrate_material=0, mask_material=2, ambient_material=3):
    """
    Generates the LevelSets and MaterialMap using the provided parameters.
    Wrapper function that creates material map and calls make_structure.

    Expects params to contain:
      - gridDelta, xExtent, yExtent (for 3D), substrateHeight, maskHeight,
        holeRadius, ambientHeight

    Returns: (list_of_levelsets, material_map, bounds, boundary_conditions, grid_delta)
    """
    dimension = int(params.get("dimensions", params.get("dimension", 3)))
    print(f"--- Generating {dimension}D Oxidation Geometry (ViennaLS) ---")

    # Set dimension for viennals
    vls.setDimension(dimension)

    # Create material map
    mat_map = vls.MaterialMap()

    # Generate structure
    level_sets = make_structure(params, mat_map, substrate_material, mask_material, ambient_material)

    # Get parameters for return values
    grid_delta = params["gridDelta"]
    x_extent = params["xExtent"]
    y_extent = params.get("yExtent", x_extent)
    substrate_height = params["substrateHeight"]
    mask_height = params["maskHeight"]
    ambient_height = params.get("ambientHeight", params.get("coverHeight", 0.0))

    z_max = substrate_height + mask_height + ambient_height + grid_delta

    if dimension == 2:
        bounds = [-x_extent/2.0, x_extent/2.0, 0.0, z_max]
        bc = [
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.INFINITE_BOUNDARY
        ]
    else:  # dimension == 3
        bounds = [-x_extent/2.0, x_extent/2.0, -y_extent/2.0, y_extent/2.0, 0.0, z_max]
        bc = [
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
            vls.BoundaryConditionEnum.INFINITE_BOUNDARY
        ]

    print(f"Geometry generated with {len(level_sets)} layers.")
    return level_sets, mat_map, bounds, bc, grid_delta

def make_cell_set(cell_set, params, substrate_material=0, mask_material=2, ambient_material=3):
    """
    Creates a DenseCellSet from level sets.
    Matches the C++ makeCellSet template function.

    Args:
        cell_set: DenseCellSet object to populate
        params: Dictionary containing geometry parameters
        substrate_material: Material ID for substrate (default: 0)
        mask_material: Material ID for mask (default: 2)
        ambient_material: Material ID for ambient (default: 3)
    """
    # Set dimension
    dimension = int(params.get("dimensions", params.get("dimension", 3)))
    vls.setDimension(dimension)

    # Generate Geometry (Level Sets)
    mat_map = vls.MaterialMap()
    level_sets = make_structure(params, mat_map, substrate_material, mask_material, ambient_material)

    # Create Cell Set (Discretization)
    depth = params.get("substrateHeight", 0.0) + params.get("ambientHeight", params.get("coverHeight", 0.0))
    cell_set.setCellSetPosition(True)  # Above surface
    cell_set.fromLevelSets(level_sets, mat_map, depth)

if __name__ == "__main__":
    import sys

    # Test with default values if run directly
    test_params = {
        "gridDelta": 0.85,
        "xExtent": 85.0,
        "yExtent": 85.0,
        "substrateHeight": 50.0,
        "maskHeight": 10.0,
        "holeRadius": 20.0,
        "ambientHeight": 8.0,
        "dimensions": 3
    }

    # Allow dimension to be specified on command line
    dimension = test_params["dimensions"]

    print(f"Testing {dimension}D geometry...")

    # Import appropriate viennacs module
    if dimension == 2:
        vls.setDimension(dimension)
        import viennacs2d as vcs
    else:
        vls.setDimension(dimension)
        import viennacs3d as vcs

    # Test 1: get_geometry_domains + manual cell set creation
    print("\nTest 1: Testing get_geometry_domains() + manual fromLevelSets...")
    lss, mat_map, _, _, _ = get_geometry_domains(test_params)

    # Also test creating a cell set from it
    cell_set1 = vcs.DenseCellSet()
    depth1 = test_params["substrateHeight"] + test_params["ambientHeight"]
    cell_set1.setCellSetPosition(True)
    print(f"  Level sets: {len(lss)}, MaterialMap layers: {mat_map.getNumberOfLayers()}")
    print(f"  Depth: {depth1}")
    cell_set1.fromLevelSets(lss, mat_map, depth1)
    cell_set1.buildNeighborhood()
    cell_set1.writeVTU(f"check_cellset_manual_{dimension}d.vtu")
    print(f"  Written check_cellset_manual_{dimension}d.vtu with {cell_set1.getNumberOfCells()} cells")

    # Also write surface mesh
    mesh = vls.Mesh()
    vls.ToSurfaceMesh(lss[-1], mesh).apply()
    vls.VTKWriter(mesh, f"check_geometry_{dimension}d.vtp").apply()
    print(f"  Written check_geometry_{dimension}d.vtp surface mesh")

    # Test 2: make_cell_set (matches C++ interface)
    print("\nTest 2: Testing make_cell_set()...")
    cell_set2 = vcs.DenseCellSet()
    make_cell_set(cell_set2, test_params,
                  substrate_material=0,
                  mask_material=2,
                  ambient_material=3)
    cell_set2.buildNeighborhood()
    cell_set2.writeVTU(f"check_cellset_{dimension}d.vtu")
    print(f"  Written check_cellset_{dimension}d.vtu with {cell_set2.getNumberOfCells()} cells")
    print("\nAll tests complete!")
