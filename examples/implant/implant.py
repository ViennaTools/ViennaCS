import sys


def read_config_file(file_name: str) -> dict:
    """Read a simple key=value config file."""
    params = {}

    with open(file_name, "r", encoding="utf-8") as file:
        for raw_line in file:
            line = raw_line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue

            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()

            try:
                params[key] = float(value)
            except ValueError:
                params[key] = value

    return params


def make_plane(x_extent: float, y_extent: float, grid_delta: float):
    try:
        import viennals as ls
        from viennals import d3 as vls
    except ImportError as exc:
        try:
            import viennals3d as vls
            ls = vls
        except ImportError:
            raise ImportError(
                "The ViennaLS Python library is required to run this example. "
                "Please install the `viennals` package and try again."
            ) from exc
    import viennacs3d as vcs

    bounds = [0.0] * (2 * vcs.D)
    bounds[1] = x_extent
    bounds[3] = y_extent

    boundary_conditions = [ls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * vcs.D
    boundary_conditions[vcs.D - 1] = ls.BoundaryConditionEnum.INFINITE_BOUNDARY

    level_set = vls.Domain(bounds, boundary_conditions, grid_delta)

    origin = [0.0] * vcs.D
    normal = [0.0] * vcs.D
    normal[vcs.D - 1] = 1.0

    vls.MakeGeometry(level_set, vls.Plane(origin, normal)).apply()
    return level_set


def main() -> int:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <config file>")
        return 1
    import viennacs3d as vcs

    params = read_config_file(sys.argv[1])

    grid_delta = params.get("gridDelta", 0.5)
    x_extent = params.get("xExtent", 10.0)
    y_extent = params.get("yExtent", 10.0)
    angle = params.get("angle", 7.0)

    level_set = make_plane(x_extent, y_extent, grid_delta)
    cell_set = vcs.DenseCellSet()
    cell_set.fromLevelSets([level_set], None, -10.0)

    cell_set.writeVTU("test.vtu")

    model = vcs.ImplantGaussian(5.0, 1.0, 2.0, 1.0)
    cell_set.addScalarData("concentration", 0.0)

    implant = vcs.Implant()
    implant.setCellSet(cell_set)
    implant.setImplantModel(model)
    implant.setImplantAngle(angle)
    implant.apply()

    cell_set.writeVTU("final.vtu")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
