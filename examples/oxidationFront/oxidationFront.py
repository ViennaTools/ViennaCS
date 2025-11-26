import os
import sys
import numpy as np
import viennacs2d as vcs
from scipy import sparse
from scipy.sparse.linalg import splu

import geometry

# ------------------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------------------
PARAMS = {
    "gridDelta": 0.85,
    "numThreads": 4,
    "substrateHeight": 50.0,
    "coverHeight": 8.0,
    "maskHeight": 10.0,
    "holeRadius": 20.0,
    "xExtent": 85.0,
    "boundaryValue": 1.0,      
    "diffusivityOxide": 2.0,
    "reactionRateConstant": 1.0,
    "eFieldInfluence": 0.5,
    "oxideConversionRate": 0.1,
    "timeStabilityFactor": 0.4,
    "duration": 300.0,
    "EfieldFile": "Efield.csv"
}

MAT_SUBSTRATE = 0
MAT_MASK = 1
MAT_COVER = 2 

class OxidationSimulation:
    def __init__(self, params):
        self.p = params
        vcs.setNumThreads(self.p['numThreads'])
        
        # 1. Geometry
        ls_list, mat_map, _, _, grid_delta = geometry.get_geometry_domains(self.p)
        self.grid_delta = grid_delta
        
        # 2. DenseCellSet
        print("Building DenseCellSet...")
        self.cell_set = vcs.DenseCellSet()
        depth = self.p['substrateHeight'] + self.p['coverHeight'] + self.p['maskHeight'] + 5.0
        
        self.cell_set.setCellSetPosition(True)
        self.cell_set.setCoverMaterial(MAT_COVER)
        self.cell_set.fromLevelSets(ls_list, mat_map, depth)
        self.cell_set.buildNeighborhood()
        self.num_cells = self.cell_set.getNumberOfCells()
        
        # 3. Data Fields
        # Register fields so C++ allocates memory for them
        self.cell_set.addScalarData("oxidant", 0.0)
        self.cell_set.addScalarData("oxideFraction", 0.0)
        # Use native Vector storage for E-field
        self.cell_set.addVectorData("Efield", [0.0, 0.0, 0.0])
        
        # Python State (Numpy Arrays for calculation)
        self.materials = np.array(self.cell_set.getScalarData("Material"), dtype=int)
        self.oxidant = np.zeros(self.num_cells)
        self.oxide_fraction = np.zeros(self.num_cells)
        self.e_field = np.zeros((self.num_cells, 3)) 
        
        # 4. Map E-Field
        self._map_electric_field()
        
        # 5. Topology
        print("Building Sparse Topology...")
        self._build_topology()
        
        # 6. Init Fields
        is_cover = self.materials == MAT_COVER
        self.oxidant[is_cover] = self.p['boundaryValue']
        
        # Sync initial state and write
        self._sync_data()
        self.cell_set.writeVTU("oxidation_initial.vtu")

    def _sync_data(self):
        """
        Copies Python numpy arrays to C++ memory. 
        Must be called before writeVTU so the file contains current data.
        """
        self.cell_set.setScalarData("oxidant", self.oxidant)
        self.cell_set.setScalarData("oxideFraction", self.oxide_fraction)
        
        # Native Vector Sync (Fast and Clean)
        self.cell_set.setVectorData("Efield", self.e_field)

    def _map_electric_field(self):
        filename = self.p['EfieldFile']
        if not os.path.exists(filename):
            print(f"[ERROR] {filename} not found.")
            sys.exit(1)
            
        try:
            # Efficient load of numeric CSV
            raw_data = np.genfromtxt(filename, delimiter=',')
        except Exception as e:
            print(f"[ERROR] Reading E-field: {e}")
            sys.exit(1)

        centers = np.array([self.cell_set.getCellCenter(i) for i in range(self.num_cells)])
        gen_x = centers[:, 0] + (self.p['xExtent'] / 2.0)
        gen_y = centers[:, 1]
        
        dx = self.grid_delta
        nx = int(self.p['xExtent'] / dx) 
        ix = np.clip((gen_x / dx).astype(int), 0, nx - 1)
        iy = (gen_y / dx).astype(int)
        
        grid_indices = iy * nx + ix
        valid_mask = (grid_indices >= 0) & (grid_indices < len(raw_data))
        
        if raw_data.shape[1] >= 3:
            self.e_field[valid_mask] = raw_data[grid_indices[valid_mask], :3]
            
        print(f"Mapped E-field to {np.sum(valid_mask)} cells.")

    def _build_topology(self):
        rows = []
        cols = []
        self.substrate_indices = np.where(self.materials == MAT_SUBSTRATE)[0]
        self.substrate_map = {idx: i for i, idx in enumerate(self.substrate_indices)}
        self.n_dof = len(self.substrate_indices)
        
        for cell_idx in self.substrate_indices:
            row_id = self.substrate_map[cell_idx]
            rows.append(row_id)
            cols.append(row_id)
            neighbors = self.cell_set.getNeighbors(cell_idx)
            for n_idx in neighbors:
                if n_idx < 0: continue
                mat_n = self.materials[n_idx]
                if mat_n == MAT_MASK: continue
                if mat_n == MAT_SUBSTRATE and n_idx in self.substrate_map:
                    col_id = self.substrate_map[n_idx]
                    rows.append(row_id)
                    cols.append(col_id)
        
        self.topo_rows = np.array(rows)
        self.topo_cols = np.array(cols)

    def _get_reaction_rate(self):
        E_mag = np.linalg.norm(self.e_field, axis=1)
        k_base = self.p['reactionRateConstant']
        alpha = self.p['eFieldInfluence']
        available_si = np.maximum(0.0, 1.0 - self.oxide_fraction)
        return k_base * (1.0 + alpha * E_mag) * available_si

    def solve_step(self, dt):
        dx = self.grid_delta
        dtdx2 = dt / (dx * dx)
        boundary_val = self.p['boundaryValue']
        D_ox = self.p['diffusivityOxide']
        
        diffusivities = np.zeros(self.num_cells)
        mask_sub = (self.materials == MAT_SUBSTRATE)
        diffusivities[mask_sub] = D_ox * self.oxide_fraction[mask_sub]
        diffusivities[self.materials == MAT_COVER] = D_ox 
        
        rhs = np.zeros(self.n_dof)
        k_rates = self._get_reaction_rate()
        
        triplet_rows = []
        triplet_cols = []
        triplet_data = []
        
        for i, cell_idx in enumerate(self.substrate_indices):
            row = i
            reaction_sink = dt * k_rates[cell_idx]
            d_sum_neighbors = 0.0
            D_center = diffusivities[cell_idx]
            val_rhs = self.oxidant[cell_idx]
            
            neighbors = self.cell_set.getNeighbors(cell_idx)
            for n_idx in neighbors:
                if n_idx < 0: continue
                mat_n = self.materials[n_idx]
                if mat_n == MAT_MASK: continue
                
                d_eff = 0.0
                if mat_n == MAT_SUBSTRATE:
                    if n_idx in self.substrate_map:
                        col = self.substrate_map[n_idx]
                        D_neighbor = diffusivities[n_idx]
                        d_eff = 0.5 * (D_center + D_neighbor)
                        triplet_rows.append(row)
                        triplet_cols.append(col)
                        triplet_data.append(-dtdx2 * d_eff)
                elif mat_n == MAT_COVER:
                    d_eff = D_center
                    val_rhs += dtdx2 * D_ox * boundary_val 
                
                d_sum_neighbors += d_eff
            
            diag_val = 1.0 + (dtdx2 * d_sum_neighbors) + reaction_sink
            triplet_rows.append(row)
            triplet_cols.append(row)
            triplet_data.append(diag_val)
            rhs[row] = val_rhs

        matrix = sparse.coo_matrix((triplet_data, (triplet_rows, triplet_cols)), 
                                   shape=(self.n_dof, self.n_dof)).tocsc()
        
        solver = splu(matrix)
        new_oxidant_sub = solver.solve(rhs)
        
        self.oxidant[self.substrate_indices] = np.clip(new_oxidant_sub, 0.0, 1.0)
        
        subset = self.substrate_indices
        C = self.oxidant[subset]
        k = k_rates[subset]
        conv_factor = self.p['oxideConversionRate']
        d_frac = k * C * dt * conv_factor
        self.oxide_fraction[subset] = np.clip(self.oxide_fraction[subset] + d_frac, 0.0, 1.0)

    def run(self):
        duration = self.p['duration']
        D_ox = self.p['diffusivityOxide']
        dt = (self.grid_delta**2 / (D_ox * 4)) * self.p['timeStabilityFactor'] * 5.0
        
        time = 0.0
        step = 0
        print(f"Starting simulation. Duration: {duration}, dt: {dt:.4f}")
        
        while time < duration:
            if time + dt > duration: dt = duration - time
            
            self.solve_step(dt)
            time += dt
            step += 1
            
            self.cell_set.writeVTU("oxidation_step_0.vtu")

            if step % 10 == 0:
                # MUST sync before writing
                self._sync_data()
                
                print(f"Step {step}, Time: {time:.4f}")
                
                self.cell_set.writeVTU(f"oxidation_step_{step}.vtu")
        
        self._sync_data()
        self.cell_set.writeVTU("oxidation_final.vtu")
        print("Simulation Done.")

if __name__ == "__main__":
    sim = OxidationSimulation(PARAMS)
    sim.run()