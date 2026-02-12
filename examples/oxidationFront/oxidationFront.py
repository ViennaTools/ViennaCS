import os
import sys
import numpy as np

import geometry

# ------------------------------------------------------------------------------
# CONFIG FILE READER - Matches ViennaPS config format
# ------------------------------------------------------------------------------
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
                par_value = line[line.find("=") + 1 :].strip()

                try:
                    val = float(par_value)
                except:
                    val = par_value

                par_dict[par_name] = val

    return par_dict

# Material IDs to match C++ version
MAT_SUBSTRATE = 0
MAT_OXIDE = 1
MAT_MASK = 2
MAT_AMBIENT = 3

class OxidationSimulation:
    """
    Electric Field Driven Oxidation Model - Python implementation.
    Matches C++ oxidationFront.cpp functionality for 2D and 3D geometry.

    Simulates Si -> SiO2 oxidation where:
    - Oxidant diffuses from ambient through oxide to silicon
    - Ambient/oxide interface has prescribed concentration (Dirichlet BC)
    - Reaction rate depends on E-field magnitude
    - Tracks oxide growth via "oxideFraction" scalar
    - Uses dynamic narrow band active cell tracking
    """
    def __init__(self, params):
        self.p = params
        self.dimension = int(params.get("dimensions", params.get("dimension", 3)))

        # Import appropriate viennacs module based on dimension
        if self.dimension == 2:
            import viennacs2d as vcs
        else:
            import viennacs3d as vcs

        vcs.setNumThreads(int(self.p['numThreads']))

        # 1. Geometry
        ls_list, mat_map, _, _, grid_delta = geometry.get_geometry_domains(self.p)
        self.grid_delta = grid_delta

        # 2. DenseCellSet
        print("Building DenseCellSet...")
        self.cell_set = vcs.DenseCellSet()
        depth = self.p['substrateHeight'] + self.p.get('ambientHeight', grid_delta)

        self.cell_set.setCellSetPosition(True)
        self.cell_set.fromLevelSets(ls_list, mat_map, depth)
        self.cell_set.buildNeighborhood()
        self.num_cells = self.cell_set.getNumberOfCells()

        # 3. Data Fields
        # Register fields so C++ allocates memory for them
        self.cell_set.addScalarData("oxidant", 0.0)
        self.cell_set.addScalarData("oxideFraction", 0.0)

        # Python State (Numpy Arrays for calculation)
        self.materials = np.array(self.cell_set.getScalarData("Material"), dtype=int)
        self.oxidant = np.zeros(self.num_cells)
        self.oxide_fraction = np.zeros(self.num_cells)
        self.e_field_1d = np.zeros(self.num_cells)  # Store scalar E-field magnitude

        # Pre-calculate neighbor map for vectorization
        self._build_neighbor_map()

        # 4. Initialize oxidant in ambient cells
        is_ambient = self.materials == MAT_AMBIENT
        self.oxidant[is_ambient] = self.p['ambientOxidant']

        # 5. Pre-compute E-Field mapping
        self._precompute_efield_indices()

        # 6. Build Topology
        print("Building Sparse Topology...")
        self._build_topology()

        # 7. Sync initial state and write
        self._sync_data()
        self.cell_set.writeVTU("oxidation_initial.vtu")

    def _build_neighbor_map(self):
        """Builds a static neighbor map for vectorized operations."""
        print("Building neighbor map for vectorization (this may take a moment)...")
        self.max_neighbors = 2 * self.dimension
        # Initialize with -1 (invalid index)
        self.neighbor_map = np.full((self.num_cells, self.max_neighbors), -1, dtype=np.int32)

        for idx in range(self.num_cells):
            nbs = self.cell_set.getNeighbors(idx)
            # Copy neighbors to array, truncating if necessary
            count = min(len(nbs), self.max_neighbors)
            if count > 0:
                self.neighbor_map[idx, :count] = nbs[:count]

    def _sync_data(self):
        """
        Copies Python numpy arrays to C++ memory.
        Must be called before writeVTU so the file contains current data.
        """
        self.cell_set.setScalarData("oxidant", self.oxidant.tolist())
        self.cell_set.setScalarData("oxideFraction", self.oxide_fraction.tolist())
        self.cell_set.setScalarData("Material", self.materials.tolist())

    def _precompute_efield_indices(self):
        """
        Pre-compute indices for E-field mapping to avoid re-calculation every step.
        """
        csv_dx = 1.0  # Should match writeCSV.py (grid of Efield source)
        csv_extent = 150.0
        nx = int(csv_extent / csv_dx)

        centers = np.array([self.cell_set.getCellCenter(i) for i in range(self.num_cells)])
        gen_x = centers[:, 0] + (csv_extent / 2.0)
        ix = np.clip((gen_x / csv_dx).astype(int), 0, nx - 1)

        # Index calculation
        if self.dimension == 2:
            # For 2D: idx = ix
            self.efield_grid_indices = ix
        else:
            # For 3D: idx = iy * nx + ix
            gen_y = centers[:, 1] + (csv_extent / 2.0)
            iy = np.clip((gen_y / csv_dx).astype(int), 0, nx - 1)
            self.efield_grid_indices = iy * nx + ix

    def _update_electric_field(self, time):
        """
        Populates self.e_field_1d with electric field magnitude values.

        This method currently reads data from a CSV file, but serves as a template for
        integrating other data sources. The process involves:
        1. Fetching E-field data (e.g., from file, database, or calculation).
        2. Mapping the external data to the simulation grid.

        The mapping uses `self.efield_grid_indices`, which is pre-computed in
        `_precompute_efield_indices` and contains the index in the external data source
        corresponding to each cell in the simulation domain.
        """
        filename = self.p['EfieldFile']
        if not os.path.exists(filename):
            print(f"[ERROR] {filename} not found.")
            sys.exit(1)

        try:
            # Load CSV - expect format matching C++ expectations
            raw_data = np.genfromtxt(filename, delimiter=',')
            # print(f"Loaded E-field data: {raw_data.shape}")
        except Exception as e:
            print(f"[ERROR] Reading E-field: {e}")
            sys.exit(1)

        # Extract E-field magnitude (matches C++ lines 171-178)
        # For D=2: expect x, E (2 cols) -> take column 1
        # For D=3: expect x, y, E (3 cols) -> take column 2
        if len(raw_data.shape) == 1:
            # Single row
            e_field_values = np.array([raw_data[-1]])
        else:
            # Multiple rows
            if self.dimension == 2 and raw_data.shape[1] >= 2:
                e_field_values = raw_data[:, 1]
            elif self.dimension == 3 and raw_data.shape[1] >= 3:
                e_field_values = raw_data[:, 2]
            else:
                # Fallback: take last column
                e_field_values = raw_data[:, -1] if raw_data.shape[1] > 1 else raw_data[:, 0]

        # Map E-field to cells
        valid_mask = (self.efield_grid_indices >= 0) & (self.efield_grid_indices < len(e_field_values))
        self.e_field_1d[:] = 0.0
        self.e_field_1d[valid_mask] = e_field_values[self.efield_grid_indices[valid_mask]]

    def _build_topology(self):
        """
        Build sparse matrix topology for active cells using narrow band approach.

        CoreActive = cells with oxideFraction > 1e-6 OR touching ambient
        Active = CoreActive OR neighbors of CoreActive (narrow band)
        """
        # Vectorized implementation
        is_mat_active = (self.materials == MAT_SUBSTRATE) | (self.materials == MAT_OXIDE)

        # Core Active Condition 1: Oxide Fraction
        cond_oxide = self.oxide_fraction > 1e-6

        # Core Active Condition 2: Touching Ambient
        valid_nb = self.neighbor_map >= 0
        safe_nbs = np.where(valid_nb, self.neighbor_map, 0)
        nb_mats = self.materials[safe_nbs]
        
        # Check if any neighbor is ambient (masking invalid neighbors)
        nb_is_ambient = (nb_mats == MAT_AMBIENT) & valid_nb
        touching_ambient = np.any(nb_is_ambient, axis=1)

        is_core_active = is_mat_active & (cond_oxide | touching_ambient)

        # Expand to Narrow Band: Active if CoreActive OR neighbor of CoreActive
        nb_is_core = is_core_active[safe_nbs] & valid_nb
        touching_core = np.any(nb_is_core, axis=1)

        is_active = is_mat_active & (is_core_active | touching_core)

        self.active_indices = np.where(is_active)[0]
        self.is_active_cell = is_active  # Store for neighbor checks in diffusion solver
        self.n_dof = len(self.active_indices)

        print(f"Active cells: {self.n_dof} / {self.num_cells} (CoreActive: {np.sum(is_core_active)})")

    def solve_step(self, dt):
        """
        Solve one time step using explicit finite difference.

        Parameters used from configuration:
        - ambientOxidant: Concentration of oxidant in the gas/ambient (Dirichlet BC).
        - oxidantDiffusivity: Diffusion coefficient of oxidant inside SiO2.
        - reactionRateConstant: Base reaction rate for Si -> SiO2 conversion.
        - eFieldInfluence: Coefficient (alpha) determining how much the E-field magnitude enhances the reaction rate.
        """
        dx = self.grid_delta
        dtdx2 = dt / (dx * dx)
        ambient_oxidant = self.p['ambientOxidant']
        D_ox = self.p['oxidantDiffusivity']

        # Vectorized Reaction Rates
        k_base = self.p['reactionRateConstant']
        alpha = self.p['eFieldInfluence']
        E_mag = np.abs(self.e_field_1d)
        available_si = np.maximum(0.0, 1.0 - self.oxide_fraction)
        reaction_rates = k_base * (1.0 + alpha * E_mag) * available_si

        # --- Diffusion (Active Cells Only) ---
        active_idxs = self.active_indices
        if len(active_idxs) > 0:
            C_old = self.oxidant[active_idxs]
            k = reaction_rates[active_idxs]
            
            # Center Diffusivity: D = D_ox * oxide_frac
            D_center = D_ox * self.oxide_fraction[active_idxs]
            
            # Neighbors
            nbs = self.neighbor_map[active_idxs] # (N_active, max_nb)
            valid_nb = nbs >= 0
            safe_nbs = np.where(valid_nb, nbs, 0)
            
            mat_nbs = self.materials[safe_nbs]
            C_nbs = self.oxidant[safe_nbs]
            ox_nbs = self.oxide_fraction[safe_nbs]
            
            diffusion_term = np.zeros_like(C_old)
            
            # Iterate over neighbor columns
            for i in range(self.max_neighbors):
                valid = valid_nb[:, i]
                mat_n = mat_nbs[:, i]
                idx_n = safe_nbs[:, i]
                
                # Ambient BC
                is_amb = (mat_n == MAT_AMBIENT) & valid
                diffusion_term[is_amb] += D_ox * (ambient_oxidant - C_old[is_amb])
                
                # Active Neighbor
                is_active_nb = self.is_active_cell[idx_n] & valid & (mat_n != MAT_AMBIENT) & (mat_n != MAT_MASK)
                
                if np.any(is_active_nb):
                    D_n = D_ox * ox_nbs[is_active_nb, i]
                    D_eff = 0.5 * (D_center[is_active_nb] + D_n)
                    diffusion_term[is_active_nb] += D_eff * (C_nbs[is_active_nb, i] - C_old[is_active_nb])

            # Explicit update
            next_oxidant_active = C_old + dtdx2 * diffusion_term - dt * k * C_old
            self.oxidant[active_idxs] = np.maximum(0.0, next_oxidant_active)

        # --- Update Oxide Fraction ---
        # Vectorized update for all substrate/oxide cells with oxidant
        mask_update = ((self.materials == MAT_SUBSTRATE) | (self.materials == MAT_OXIDE)) & (self.oxidant > 1e-12)
        idxs_update = np.where(mask_update)[0]
        
        if len(idxs_update) > 0:
            C = self.oxidant[idxs_update]
            k = reaction_rates[idxs_update]
            d_frac = k * C * dt
            
            self.oxide_fraction[idxs_update] += d_frac
            self.oxide_fraction[idxs_update] = np.minimum(1.0, self.oxide_fraction[idxs_update])
            
            # Update Material if > 50% oxidized
            new_oxides = (self.oxide_fraction[idxs_update] > 0.5) & (self.materials[idxs_update] == MAT_SUBSTRATE)
            if np.any(new_oxides):
                self.materials[idxs_update[new_oxides]] = MAT_OXIDE

    def run(self):
        """
        Main simulation loop with adaptive time stepping.
        """
        duration = self.p['duration']
        D_ox = self.p['oxidantDiffusivity']

        # Estimate stable time step
        # For 2D: 2*D = 4, For 3D: 2*D = 6
        stability_factor = 2 * self.dimension
        dt_base = (self.grid_delta**2 / (D_ox * stability_factor)) * self.p['timeStabilityFactor']

        time = 0.0
        step = 0
        print(f"Starting simulation. Duration: {duration}, initial dt: {dt_base:.4f}")

        while time < duration:
            # Update active cells
            self._build_topology()

            # Update E-field
            self._update_electric_field(time)

            # Adaptive time stepping based on max reaction rate
            # dt depends on:
            # 1. E-field at active cells (spatially varying)
            # 2. Available silicon (decreases as oxidation proceeds)
            if len(self.active_indices) > 0:
                # Calculate k for active cells to determine dt
                k_base = self.p['reactionRateConstant']
                alpha = self.p['eFieldInfluence']
                E_mag = np.abs(self.e_field_1d[self.active_indices])
                avail_si = np.maximum(0.0, 1.0 - self.oxide_fraction[self.active_indices])
                rates = k_base * (1.0 + alpha * E_mag) * avail_si
                max_k = np.max(rates)
            else:
                max_k = 0.0

            dt_diff = self.grid_delta**2 / (D_ox * stability_factor)
            dt_react = (1.0 / max_k) if max_k > 1e-12 else duration
            dt = min(dt_diff, dt_react) * self.p['timeStabilityFactor']

            if time + dt > duration:
                dt = duration - time

            self.solve_step(dt)
            time += dt
            step += 1

            if step % 10 == 0:
                # Sync and write output
                self._sync_data()
                print(f"Step {step}, Time: {time:.4f}, dt: {dt:.4f}, max_k: {max_k:.4f}")
                self.cell_set.writeVTU(f"oxidation_step_{step}.vtu")

        self._sync_data()
        self.cell_set.writeVTU("oxidation_final.vtu")
        print("Simulation Done.")

if __name__ == "__main__":
    # Read config file from command line argument
    # Usage: oxidationFront.py [config_file]
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <config file>")
        print("Using default config.txt")
        config_file = "config.txt"
    else:
        config_file = sys.argv[1]

    if not os.path.exists(config_file):
        print(f"[ERROR] Config file '{config_file}' not found.")
        sys.exit(1)

    params = ReadConfigFile(config_file)

    # Get dimension from config file
    dimension = int(params.get("dimensions", params.get("dimension", 3)))

    print(f"Loaded configuration from {config_file}")
    print(f"  Dimension: {dimension}D")
    print(f"  Grid delta: {params['gridDelta']}")
    if dimension == 2:
        print(f"  Domain: {params['xExtent']} x Y")
    else:
        y_extent = params.get('yExtent', params['xExtent'])
        print(f"  Domain: {params['xExtent']} x {y_extent} x Z")
    print(f"  Duration: {params['duration']}")

    sim = OxidationSimulation(params)
    sim.run()
