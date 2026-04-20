import viennacs as vcs
import viennals as vls

D = 2
vcs.setDimension(D)
vls.setDimension(D)


# Config file reader helper function
def ReadConfigFile(fileName: str):
    """Read a config file in the ViennaPS standard config file format.

    Parameters
    ----------
    fileName: str
              Name of the config file.

    Returns
    -------
    dict
        A dictionary containing the parameters from the config file.
    """
    par_dict = {}

    with open(fileName, "r") as file:
        lines = file.readlines()
        for line in lines:

            line = line[: line.find("#")]  # remove comments

            if len(line) > 0:
                par_name = line[: line.find("=")].strip(" ")
                par_value = line[line.find("=") + 1 :]

                try:
                    val = float(par_value)
                except:
                    val = par_value

                par_dict[par_name] = val

    return par_dict


def makeLShape(params: dict):
    try:
        import viennals as vls
    except ImportError:
        print(
            "The ViennaLS Python library is required to run this example. "
            "Please install it and try again."
        )
        return

    gridDelta = params["gridDelta"]
    bounds = [0] * D * 2
    bounds[0] = -params["verticalWidth"] / 2.0 - params["xPad"]
    bounds[1] = params["verticalWidth"] / 2.0 + params["xPad"]
    if D == 3:
        bounds[2] = -params["verticalWidth"] / 2.0 - params["xPad"]
        bounds[3] = (
            -params["verticalWidth"] / 2.0 + params["xPad"] + params["horizontalWidth"]
        )
    else:
        bounds[1] = (
            -params["verticalWidth"] / 2.0 + params["xPad"] + params["horizontalWidth"]
        )

    boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * (D - 1)
    boundaryCons.append(vls.BoundaryConditionEnum.INFINITE_BOUNDARY)

    substrate = vls.Domain(bounds, boundaryCons, gridDelta)
    normal = [0.0] * D
    origin = [0.0] * D
    normal[D - 1] = 1.0
    origin[D - 1] = params["verticalDepth"]
    vls.MakeGeometry(substrate, vls.Plane(origin, normal)).apply()

    # Create the vertical trench
    vertBox = vls.Domain(substrate)
    minPoint = [0] * D
    maxPoint = [0] * D
    for i in range(D - 1):
        minPoint[i] = -params["verticalWidth"] / 2.0
        maxPoint[i] = params["verticalWidth"] / 2.0
    maxPoint[D - 1] = params["verticalDepth"]
    vls.MakeGeometry(vertBox, vls.Box(minPoint, maxPoint)).apply()
    vls.BooleanOperation(
        substrate, vertBox, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
    )

    # Create the horizontal trench
    horiBox = vls.Domain(substrate)
    minPoint = [0] * D
    maxPoint = [params["verticalWidth"] / 2.0] * D
    for i in range(D - 1):
        minPoint[i] = -params["verticalWidth"] / 2.0
    maxPoint[D - 1] = params["horizontalHeight"]
    maxPoint[D - 2] = -params["verticalWidth"] / 2.0 + params["horizontalWidth"]
    vls.MakeGeometry(horiBox, vls.Box(minPoint, maxPoint)).apply()
    vls.BooleanOperation(
        substrate, horiBox, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
    )

    return substrate


# vcs.Logger.setLogLevel(vcs.LogLevel.INTERMEDIATE)
params = ReadConfigFile("config.txt")
# vps.setNumThreads(int(params["numThreads"]))

# Create a cellSet
ls = makeLShape(params)
cellSet = vcs.DenseCellSet()
cellSet.setCellSetPosition(True)
cellSet.setCoverMaterial(1)
matMap = vls.MaterialMap()
matMap.insertNextMaterial(0)
cellSet.fromLevelSets([ls], matMap, params["topSpace"] + params["verticalDepth"])

# Segment the cells into surface, material, and gas cells
vcs.SegmentCells(cellSet, "CellType", 1).apply()

# Calculate the mean free path for the gas cells
mfpCalc = vcs.MeanFreePath(cellSet)
mfpCalc.setNumRaysPerCell(params["raysPerCell"])
mfpCalc.setReflectionLimit(int(params["reflectionLimit"]))
mfpCalc.setRngSeed(int(params["seed"]))
mfpCalc.setMaterial(1)
mfpCalc.setBulkLambda(params["bulkLambda"])
mfpCalc.apply()

cellSet.writeVTU("initial.vtu")

# Run the atomic layer deposition model
model = vcs.AtomicLayerProcess(cellSet)
model.setMaxLambda(params["bulkLambda"])
model.setPrintInterval(params["printInterval"])
model.setStabilityFactor(params["stabilityFactor"])
model.setFirstPrecursor(
    "H2O",
    params["H2O_meanThermalVelocity"],
    params["H2O_adsorptionRate"],
    params["H2O_desorptionRate"],
    params["p1_time"],
    params["inFlux"],
)
model.setSecondPrecursor(
    "TMA",
    params["TMA_meanThermalVelocity"],
    params["TMA_adsorptionRate"],
    params["TMA_desorptionRate"],
    params["p2_time"],
    params["inFlux"],
)
model.setPurgeParameters(params["purge_meanThermalVelocity"], params["purge_time"])
model.setReactionOrder(params["reactionOrder"])
model.apply()
