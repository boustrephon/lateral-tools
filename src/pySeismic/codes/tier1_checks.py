"""Tier 1 capacity functions per ASCE 41-17, Section 4.4.3 and checklists."""

# from typing import Union
import math

from sys import path as sys_path
from os.path import abspath, join as path_join, dirname

from typing import Optional

# Add the parent directory of 'pySeismic' to sys.path
sys_path.append(abspath(path_join(dirname(__file__), '..', '..')))

from pySeismic.codes.ASCE7_16 import site_coeffs, calculate_k_factor, seismic_design_category_get

# Unit conversion functions (defined in ASCE41_17.py)
from pySeismic.codes.ASCE41_17 import kN2kip, sqm2sqin, MPa2psi, m2ft

# Demand functions used in wrapper checks
from pySeismic.codes.ASCE41_17 import seismicity_level_get, v_j_avg_rcmf, v_j_avg_sw, p_ot

from pySeismic.codes.ASCE41_17 import (
    table_3_4_get, # Limitations
    T_code_get, # Code period in seconds
    S_a_Tier_1_func, # Spectral acceleration
    table_4_7_dict_17_get, # C-factor
)

from pySeismic.utilities import storey_force_distribution, storey_shear_moment



# Optional: if you need the M_s table for URM infill
# from pySeismic.codes.ASCE41_17 import table_4_8_dict_17_get

# ----------------------------------------------------------------------
# Concrete material capacities (imperial units: psi, ksi)
# ----------------------------------------------------------------------

def capacity_shear_stress_rcmf(perf: str, f_c: float) -> float:
    """
    Allowable shear stress in RC moment frame columns (psi).
    
    From ASCE 41-17 Tier 1 checklists:
    - IO: 2√(f'c) [psi], minimum 100 psi (Table 17-23)
    - LS: 3√(f'c) [psi] (derived from similar provisions; adjust as needed)
    - CP: 4√(f'c) [psi] (derived from similar provisions)
    
    Args:
        perf: Performance level ('IO', 'LS', 'CP')
        f_c: Concrete compressive strength (psi)
    
    Returns:
        Allowable shear stress (psi)
    """
    f_c_psi = f_c
    if perf.upper() == 'IO':
        return max(2.0 * math.sqrt(f_c_psi), 100.0)
    elif perf.upper() == 'LS':
        return 3.0 * math.sqrt(f_c_psi)
    elif perf.upper() == 'CP':
        return 4.0 * math.sqrt(f_c_psi)
    else:
        raise ValueError(f"Unknown performance level: {perf}")


def capacity_column_axial_stress(perf: str, f_c: float) -> float:
    """
    Allowable axial stress in columns due to overturning (psi).
    
    From ASCE 41-17 Tier 1 checklists:
    - IO: 0.3 f'c (Table 17-23)
    - LS: 0.5 f'c (derived)
    - CP: 0.7 f'c (derived)
    
    Args:
        perf: Performance level ('IO', 'LS', 'CP')
        f_c: Concrete compressive strength (psi)
    
    Returns:
        Allowable axial stress (psi)
    """
    f_c_psi = f_c
    if perf.upper() == 'IO':
        return 0.3 * f_c_psi
    elif perf.upper() == 'LS':
        return 0.5 * f_c_psi
    elif perf.upper() == 'CP':
        return 0.7 * f_c_psi
    else:
        raise ValueError(f"Unknown performance level: {perf}")


def capacity_shear_stress_wall(perf: str, f_c: float) -> float:
    """
    Allowable shear stress in RC shear walls (psi).
    
    From ASCE 41-17 Tier 1 checklists (Table 17-25):
    - IO: 2√(f'c) [psi], minimum 100 psi
    - LS: 3√(f'c) [psi]
    - CP: 4√(f'c) [psi]
    
    Args:
        perf: Performance level ('IO', 'LS', 'CP')
        f_c: Concrete compressive strength (psi)
    
    Returns:
        Allowable shear stress (psi)
    """
    return capacity_shear_stress_rcmf(perf, f_c)  # same formula as columns


def capacity_joint_shear_stress(perf: str, f_c: float) -> float:
    """
    Allowable joint shear stress for beam‑column joints (psi).
    
    Per ACI 318 and ASCE 41-17 Tier 1 (e.g., Table 17-23):
    - IO: 12√(f'c) [psi] (typical for confined joints)
    - LS: 15√(f'c) [psi]
    - CP: 20√(f'c) [psi]
    
    Args:
        perf: Performance level ('IO', 'LS', 'CP')
        f_c: Concrete compressive strength (psi)
    
    Returns:
        Allowable joint shear stress (psi)
    """
    f_c_psi = f_c
    if perf.upper() == 'IO':
        return 12.0 * math.sqrt(f_c_psi)
    elif perf.upper() == 'LS':
        return 15.0 * math.sqrt(f_c_psi)
    elif perf.upper() == 'CP':
        return 20.0 * math.sqrt(f_c_psi)
    else:
        raise ValueError(f"Unknown performance level: {perf}")


def check_shear_stress_rcmf(
    V_j: float,           # story shear (kip)
    A_c: float,           # total column area (in²)
    n_c: int,             # number of columns
    n_f: int,             # number of frames in direction
    perf: str,            # performance level
    f_c: float,           # concrete strength (psi)
    M_s: float|None = None     # optional Ms; if None, uses default from ASCE41_17.v_j_avg_rcmf
) -> dict:
    """
    Perform Tier 1 shear stress check for RC moment frames.
    
    Returns a dictionary with demand, capacity, D/C ratio, and pass/fail.
    """
    # from pySeismic.codes.ASCE41_17 import v_j_avg_rcmf
    demand_psi = v_j_avg_rcmf(V_j, A_c, n_c, n_f, perf)
    capacity_psi = capacity_shear_stress_rcmf(perf, f_c)
    dc = demand_psi / capacity_psi
    return {
        'check': 'Shear Stress (columns)',
        'demand_psi': demand_psi,
        'capacity_psi': capacity_psi,
        'dc_ratio': dc,
        'pass': dc <= 1.0
    }


def check_column_axial_overturning(
    V: float,             # base shear (kip)
    h_n: float,           # building height (ft)
    L: float,             # frame length (ft)
    n_f: int,             # number of frames
    A_col: float,         # area of end column (in²)
    perf: str,
    f_c: float,
) -> dict:
    # from pySeismic.codes.ASCE41_17 import p_ot
    demand_psi = p_ot(V, h_n, L, n_f, A_col, perf)
    capacity_psi = capacity_column_axial_stress(perf, f_c)
    dc = demand_psi / capacity_psi
    return {
        'check': 'Column Axial Stress (overturning)',
        'demand_psi': demand_psi,
        'capacity_psi': capacity_psi,
        'dc_ratio': dc,
        'pass': dc <= 1.0
    }


def check_shear_stress_wall(
    V_j: float,
    A_w: float,           # total wall area in direction (in²)
    wtype: str,           # wall type ('RC' or 'Masonry')
    perf: str,
    f_c: float,
) -> dict:
    # from pySeismic.codes.ASCE41_17 import v_j_avg_sw
    demand_psi = v_j_avg_sw(V_j, A_w, wtype, perf)
    capacity_psi = capacity_shear_stress_wall(perf, f_c)
    dc = demand_psi / capacity_psi
    return {
        'check': 'Shear Stress (walls)',
        'demand_psi': demand_psi,
        'capacity_psi': capacity_psi,
        'dc_ratio': dc,
        'pass': dc <= 1.0
    }


class Tier1Checker:
    """
    Run Tier 1 screening checks for a given building.
    
    The building data should be a dictionary loaded from JSON with keys:
        'structure_type': 'C1', 'C2', or 'C3'
        'performance_level': 'IO', 'LS', 'CP'
        'seismicity_level': 'Very Low', 'Low', 'Moderate', 'High'
        'directions': {'X': {...}, 'Y': {...}} containing:
            'story_shear_kN': list of story shears per level (top to bottom)
            'total_column_area_m2': total column cross‑section area
            'num_columns': number of columns
            'num_frames': number of frames in this direction
            'frame_length_m': total frame length
            'end_column_area_m2': area of end column
            'shear_wall_area_m2': total wall area (if C2)
            'f_c_MPa': concrete strength
        'storey_elevations_m': list of elevations (top to bottom)
        'storey_weights_kN': list of weights per storey (top to bottom)
    """
    
    def __init__(self, building_data: dict, calculate_forces: bool = True, debug:bool=True):
        self.data = building_data
        self.struct_type = building_data['structure_type']
        self.perf = building_data['performance_level']
        Ss = self.data.get('Ss')
        S1 = self.data.get('S1')
        site_class = self.data.get('site_class')
        risk_category = self.data.get('risk_category')
        self.debug = debug
        if debug:
            print(f'Ss: {Ss}, S1: {S1}, site_class: {site_class}, risk_category: {risk_category}')

        # Determine S_DS and S_D1
        self._set_S_DS_S_D1()

        # Determine seismic design category (it independently calculates S_DS & S_D1 if they have not already been calculated)
        self.sdc = seismic_design_category_get(Ss, S1, site_class, risk_category, debug=True)

        # Determine seismicity level (adds to data if not present)
        self._determine_seismicity()
        
        self.seis = building_data['seismicity_level']
        print(f'Seismicity: {self.seis}')
        print(f'S_DS = {self.data["S_DS"]}, S_D1 = {self.data["S_D1"]} (init)')
    
        # Calculate seismic forces if requested
        if calculate_forces:
            self.calculate_seismic_forces()
    
        self.results = []   # list of check results


    def calculate_seismic_forces(self):
        """
        Calculate building period, base shear, and lateral force distribution.
        Updates self.data with results.
        """
        # from pySeismic.codes.ASCE41_17 import (
        #     building_period_approx_si,  # We'll use this instead of direct calculation
        #     base_shear_calculation,
        #     get_R_value,
        #     m2ft,
        #     calculate_k_factor
        # )
        # from pySeismic.utilities import storey_force_distribution, storey_shear_moment
        
        # Get required data
        struct_type = self.struct_type
        building_height_m = self.data.get('building_height_m')
        if not building_height_m: #  and ('storey_elevations_m' in self.data):
            # Calculate height from top storey elevation
            building_height_m = max(self.data['storey_elevations_m'])
            self.data['building_height_m'] = building_height_m
        num_storeys = self.data.get('num_storeys')
        if not num_storeys: #  and ('num_storeys' in self.data):
            # Calculate storeys from storey weights
            num_storeys = len(self.data['storey_elevations_m'])
            self.data['num_storeys'] = num_storeys
        if self.debug:
            print(f'struct_type: {struct_type}, building_height_m: {building_height_m}, num_storeys: {num_storeys}')
        
        # Get spectral values (already computed in _determine_seismicity)
        S_DS = self.data.get('S_DS')
        S_D1 = self.data.get('S_D1')
        if self.debug:
            print(f'S_DS: {S_DS}, S_D1: {S_D1} (repeat)')
        
        if not all([building_height_m, num_storeys, S_DS, S_D1]):
            missing = []
            if not building_height_m: 
                missing.append('building_height_m')
            if not num_storeys: 
                missing.append('num_storeys')
            if not S_DS: 
                missing.append('S_DS')
            if not S_D1: 
                missing.append('S_D1')
            raise ValueError(f"Cannot calculate seismic forces: missing {missing}")
        
        # Calculate period
        T_code = T_code_get(struct_type, building_height_m, 'm')
        C = table_4_7_dict_17_get(struct_type, num_storeys)
        S_a = S_a_Tier_1_func(T_code, self.data['Ss'], self.data['S1'], self.data['site_class'])
        if self.debug:
            print(f'T_code: {T_code:.2f} sec, C: {C:.2f}, S_a: {S_a:.3f}g')
        Cs = C * S_a
        
        # Store period results
        self.data['period_analysis'] = {
            'T_code': T_code,
            'C': C,
            'S_a': S_a,
            'Cs': Cs
        }
        self.data['C'] = C
        
        # Calculate total building weight
        storey_weights_kN = self.data.get('storey_weights_kN', [])
        if not storey_weights_kN:
            raise ValueError("storey_weights_kN required for base shear calculation")
        
        total_weight_kN = sum(storey_weights_kN)
        self.data['total_weight_kN'] = total_weight_kN
        
        # Calculate base shear in kN
        V_base_kN = Cs * total_weight_kN
        self.data['base_shear_kN'] = V_base_kN
        
        # Calculate lateral force distribution
        storey_elevations_m = self.data.get('storey_elevations_m', [])
        if not storey_elevations_m:
            raise ValueError("storey_elevations_m required for force distribution")
        
        # Ensure elevations are top-to-bottom (descending) for storage,
        # but force distribution needs ascending order
        if storey_elevations_m[0] < storey_elevations_m[-1]:
            # Already ascending - reverse for storage
            storey_elevations_asc = storey_elevations_m
            storey_weights_asc = storey_weights_kN
            storey_elevations_desc = storey_elevations_m[::-1]
            storey_weights_desc = storey_weights_kN[::-1]
        else:
            # Descending - reverse for calculation
            storey_elevations_asc = storey_elevations_m[::-1]
            storey_weights_asc = storey_weights_kN[::-1]
            storey_elevations_desc = storey_elevations_m
            storey_weights_desc = storey_weights_kN
        
        # Calculate k factor
        k = calculate_k_factor(T_code)
        
        # Calculate force distribution
        storey_forces_asc = storey_force_distribution(
            V_base_kN, 
            storey_elevations_asc, 
            storey_weights_asc, 
            k
        )
        
        # Calculate storey shears and moments
        # First storey height is from base to first storey
        storey_heights_asc = [storey_elevations_asc[0]] + [
            storey_elevations_asc[i] - storey_elevations_asc[i-1] 
            for i in range(1, len(storey_elevations_asc))
        ]
        
        storey_shears_asc, storey_moments_asc = storey_shear_moment(
            storey_forces_asc, 
            storey_heights_asc, 
            ascending=True
        )
        
        # Convert back to top-to-bottom order for storage
        self.data['storey_forces_kN'] = storey_forces_asc[::-1]
        self.data['storey_shears_kN'] = storey_shears_asc[::-1]
        self.data['storey_moments_kNm'] = storey_moments_asc[::-1]
        self.data['k_factor'] = k
        
        # Calculate for both directions if not symmetric
        for dirn in ('X', 'Y'):
            if dirn in self.data.get('directions', {}):
                self.data['directions'][dirn]['story_shear_kN'] = self.data['storey_shears_kN']
        
        return self.data['period_analysis']


    def run_checks(self):
        """Run all checks applicable to the structure type."""
        # Basic configuration checks (example – can be expanded)
        self._check_basic_config()
        
        # Run the appropriate checklist items
        self._run_checklist_items()
    
        # Structural checks per type
        if self.struct_type in ('C1', 'C2', 'C3'):
            self._check_shear_stress_columns()
            self._check_column_axial_overturning()
        if self.struct_type == 'C2':
            self._check_shear_stress_walls()
        if self.struct_type == 'C1':
            self._check_joint_shear_stress()   # to be implemented
        if self.struct_type == 'C3':
            self._check_infill_panel_shear()   # to be implemented

        # C2-specific checks
        if self.struct_type in ('C2', 'C2a'):
            # Wall shear stress check (already implemented)
            self._check_shear_stress_walls()
            
            # New C2 methods
            self._check_wall_redundancy()
            self._check_wall_reinforcing()
            self._check_wall_aspect_ratios()
            self._check_boundary_confinement()
            self._check_wall_openings()
            self._check_wall_thickness()
            self._check_coupling_beams()
            self._check_secondary_frames()
            
            # Diaphragm checks
            self._check_diaphragm_continuity()
            self._check_diaphragm_openings_at_walls()
            self._check_diaphragm_cross_ties()
            self._check_diaphragm_sheathing()
            self._check_metal_deck()
            
            # Foundation checks
            self._check_foundation_dowels()
            self._check_pile_cap_uplift()
        
        return self.results       

    def _run_checklist_items(self):
        """
        Run all applicable checklist items for the building type and performance level.
        """
        struct_type = self.struct_type
        perf = self.perf
        seismicity = self.seis
        
        # Get the appropriate checklist
        checklist = TIER1_CHECKLISTS.get(struct_type, {}).get(perf, {})
        if not checklist:
            self.results.append({
                'check': f'No checklist found for {struct_type} {perf}',
                'pass': None,
                'note': 'Checklist not available'
            })
            return
        
        # Run each checklist item
        for item_id, item_data in checklist.items():
            # Check if this item applies to the current seismicity level
            if 'seismicity' in item_data:
                if seismicity not in item_data['seismicity']:
                    continue
            
            # Check if this item applies to the diaphragm type (if specified)
            if 'diaphragm_type' in item_data:
                # You would need to know the diaphragm type from building data
                # This is just a placeholder - you'd need to add diaphragm_type to building_data
                pass
            
            # For now, just add the statement as a check that needs manual review
            # In a real implementation, you'd have logic to evaluate each statement
            self.results.append({
                'check': item_data['statement'],
                'pass': None,  # Would need actual evaluation logic
                'tier2_ref': item_data.get('tier2_ref', ''),
                'note': 'Manual review required'
            })

    def _set_S_DS_S_D1(self):
        """Set S_DS and S_D1 from Ss and S1."""
        # Check for S_DS and S_D1 directly
        if 'S_DS' in self.data and 'S_D1' in self.data:
            if self.debug:
                print("A) 'S_DS' in self.data and 'S_D1' in self.data")
            s_ds = self.data['S_DS']
            s_d1 = self.data['S_D1']
            if self.debug:
                print(f'A) S_DS = {s_ds} | {self.data["S_DS"]}, S_D1 = {s_d1} | {self.data["S_D1"]}')

        else:
            # Need Ss, S1, and site_class
            if self.debug:
                print("B) 'S_DS' NOT in self.data or 'S_D1' NOT in self.data")
                print(f'B) Ss = {self.data["Ss"]}, S1 = {self.data["S1"]}, site_class = {self.data["site_class"]}')
            required = ('Ss', 'S1', 'site_class')
            if not all(k in self.data for k in required):
                missing = [k for k in required if k not in self.data]
                raise ValueError(f"Cannot determine seismicity: missing {missing}")
            
            # from pySeismic.codes.ASCE7_16 import site_coeffs
            # from pySeismic.codes.ASCE41_17 import seismicity_level_get

            Fa, Fv = site_coeffs(self.data['site_class'], self.data['Ss'], self.data['S1'])
            s_ds = Fa * self.data['Ss']
            s_d1 = Fv * self.data['S1']
            # Optionally store them for later use
            self.data['S_DS'] = s_ds
            self.data['S_D1'] = s_d1
            if self.debug:
                print(f'C) S_DS = {s_ds} | {self.data["S_DS"]}, S_D1 = {s_d1} | {self.data["S_D1"]}')


    def _determine_seismicity(self):
        """
        Determine seismicity level from Ss and S1, or from S_DS/S_D1 if directly provided.
        Updates self.data['seismicity_level'].
        """
        if self.debug:
            print("Determine seismicity level")
            
        # If already present, use it
        if 'seismicity_level' in self.data:
            return

        self._set_S_DS_S_D1(self, debug=False)
        s_ds = self.data['S_DS']
        s_d1 = self.data['S_D1']
        if debug:
            print(f'D) S_DS = {s_ds} | {self.data["S_DS"]}, S_D1 = {s_d1} | {self.data["S_D1"]}')

        # from pySeismic.codes.ASCE41_17 import seismicity_level_get
        seismicity = seismicity_level_get(s_ds, s_d1)
        self.data['seismicity_level'] = seismicity

    def _check_shear_stress_columns(self):
        """Run shear stress check in both orthogonal directions."""
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            # Convert SI to imperial
            V_j = kN2kip(d['story_shear_kN'][-1])  # use base shear
            A_c = sqm2sqin(d['total_column_area_m2'])
            n_c = d['num_columns']
            n_f = d['num_frames']
            f_c_psi = MPa2psi(d['f_c_MPa'])
            res = check_shear_stress_rcmf(V_j, A_c, n_c, n_f, self.perf, f_c_psi)
            res['direction'] = dirn
            self.results.append(res)
    
    def _check_column_axial_overturning(self):
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            V = kN2kip(self.data['base_shear_kN'])   # need base shear stored
            h_n = m2ft(self.data['building_height_m'])
            L = m2ft(d['frame_length_m'])
            n_f = d['num_frames']
            A_col = sqm2sqin(d['end_column_area_m2'])
            f_c_psi = MPa2psi(d['f_c_MPa'])
            res = check_column_axial_overturning(V, h_n, L, n_f, A_col, self.perf, f_c_psi)
            res['direction'] = dirn
            self.results.append(res)
    
    def _check_shear_stress_walls(self):
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d or 'shear_wall_area_m2' not in d:
                continue
            V_j = kN2kip(d['story_shear_kN'][-1])
            A_w = sqm2sqin(d['shear_wall_area_m2'])
            f_c_psi = MPa2psi(d['f_c_MPa'])
            res = check_shear_stress_wall(V_j, A_w, 'RC', self.perf, f_c_psi)
            res['direction'] = dirn
            self.results.append(res)
    
    # ----------------------------------------------------------------------
    # Basic configuration checks (Section 17.1.2)
    # ----------------------------------------------------------------------
    def _check_basic_config(self):
        """
        Evaluate configuration‑related statements from the Basic Configuration Checklist.
        The building data should contain a 'config_checks' dictionary with keys matching
        the statement IDs below. If a key is missing, the check is marked as "Not Checked".
        """
        # List of (check_id, description, default_pass_if_missing)
        config_statements = [
            ('load_path',          'Load path: A complete lateral load path exists.'),
            ('adjacent_buildings', 'Adjacent buildings: Sufficient separation exists.'),
            ('meandering',         'Meandering: Columns are aligned vertically.'),
            ('vertical_irregular', 'Vertical irregularities: No soft or weak storey.'),
            ('plan_irregular',     'Plan irregularities: No torsional irregularity.'),
            ('diaphragm_continuity','Diaphragm continuity: Diaphragms are continuous.'),
        ]
        cfg = self.data.get('config_checks', {})
        for cid, desc in config_statements:
            # If the data provides a boolean, use it; otherwise mark as None (not checked)
            passed = cfg.get(cid)
            if passed is None:
                result = {
                    'check': desc,
                    'pass': None,
                    'note': 'Data not provided – manual review required'
                }
            else:
                result = {
                    'check': desc,
                    'pass': passed,
                    'note': ''
                }
            self.results.append(result)

    # ----------------------------------------------------------------------
    # Joint shear stress for C1 (RC moment frames)
    # ----------------------------------------------------------------------
    def _check_joint_shear_stress(self):
        """
        Beam‑column joint shear stress check (C1 only).
        Requires additional data under 'joints' in each direction, containing:
            - 'joint_width_m': effective joint width (m)
            - 'beam_depth_m': average beam depth (m)   (optional; if missing uses default)
            - 'column_depth_m': column depth in direction of loading (m)
        If not provided, the check is skipped with a warning.
        """
        if self.struct_type != 'C1':
            return

        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            joints_data = d.get('joints')
            if not joints_data:
                self.results.append({
                    'check': f'Joint Shear Stress ({dirn})',
                    'pass': None,
                    'note': 'Joint data missing – cannot evaluate'
                })
                continue

            # Demand: simplified joint shear based on column shear
            # Use the largest storey shear (usually base) for a conservative estimate
            V_story_kN = max(d['story_shear_kN'])          # kN
            n_col = d['num_columns']
            V_col_kN = V_story_kN / n_col                  # average column shear

            # Convert to imperial
            V_col_kip = kN2kip(V_col_kN)
            joint_width_in = m2ft(joints_data['joint_width_m']) * 12   # convert m to inches
            column_depth_in = m2ft(joints_data.get('column_depth_m', 0.5)) * 12
            # Assume beam depth ~ 0.6 * column depth if not given
            beam_depth_in = m2ft(joints_data.get('beam_depth_m', 0.6 * joints_data.get('column_depth_m', 0.5))) * 12

            # Simplified joint shear force (kip) – based on beam tension/compression
            # Assuming M_b = V_col * (story_height/2) and jd ≈ 0.9*d_b
            story_height_ft = m2ft(self.data['storey_elevations_m'][0] - self.data['storey_elevations_m'][1])  # rough
            M_b_kft = V_col_kip * (story_height_ft / 2)                    # beam moment
            T_kip = M_b_kft * 12 / (0.9 * beam_depth_in)                   # tension force
            V_j_kip = T_kip * (column_depth_in / beam_depth_in)            # very approximate

            # Joint shear stress (psi)
            demand_psi = V_j_kip * 1000 / (joint_width_in * column_depth_in)

            # Capacity
            f_c_psi = MPa2psi(d['f_c_MPa'])
            capacity_psi = capacity_joint_shear_stress(self.perf, f_c_psi)

            dc = demand_psi / capacity_psi
            self.results.append({
                'check': f'Joint Shear Stress ({dirn})',
                'demand_psi': demand_psi,
                'capacity_psi': capacity_psi,
                'dc_ratio': dc,
                'pass': dc <= 1.0,
                'note': 'Simplified Tier 1 approximation – use with caution'
            })

    # ----------------------------------------------------------------------
    # Infill panel shear for C3 (RC frames with URM infill)
    # ----------------------------------------------------------------------
    def _check_infill_panel_shear(self):
        """
        Shear stress in URM infill panels (C3 only).
        Uses wall information from building data, filtered by orientation and type='URM'.
        For each storey, total infill wall area is summed, and average shear stress
        is compared to allowable URM shear stress.
        """
        if self.struct_type != 'C3':
            return

        # Need material properties for masonry
        mat = self.data.get('material_properties', {})
        f_m_MPa = mat.get('f_m_MPa', 2.0)          # default 2 MPa if missing
        f_m_psi = MPa2psi(f_m_MPa)

        # Group walls by storey and orientation
        walls = self.data.get('wall_information', [])
        # Assume each wall dict has: 'storey' (int, 0=top), 'orientation', 'thickness_m', 'length_m', 'type'
        from collections import defaultdict
        area_by_storey_dir = defaultdict(lambda: defaultdict(float))
        for w in walls:
            if w.get('type', '').upper() != 'URM':
                continue
            storey = w['storey']
            orient = w['orientation'].upper()
            if orient in ('X', 'Y'):
                area = w['thickness_m'] * w['length_m']   # cross‑sectional area in m²
                area_by_storey_dir[storey][orient] += area

        # For each direction, compute demand per storey (using storey shear at that level)
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            story_shears_kN = d['story_shear_kN']          # list, top to bottom
            n_storeys = len(story_shears_kN)

            for storey_idx in range(n_storeys):
                V_kN = story_shears_kN[storey_idx]
                A_w_m2 = area_by_storey_dir[storey_idx].get(dirn, 0.0)
                if A_w_m2 == 0:
                    # No infill in this storey/direction – skip or note?
                    continue

                # Convert to imperial
                V_kip = kN2kip(V_kN)
                A_w_in2 = sqm2sqin(A_w_m2)

                demand_psi = V_kip * 1000 / A_w_in2        # average shear stress

                # Allowable shear stress for URM infill (simplified)
                # ASCE 41-17 Table 4-8 gives M_s factors, but base shear strength is needed.
                # For Tier 1, we use a fraction of masonry strength. Common: 0.1*f'm to 0.2*f'm.
                capacity_psi = 0.15 * f_m_psi              # rough estimate – adjust as needed
                # Optionally use M_s factor from Table 4-8 (requires performance level)
                # from pySeismic.codes.ASCE41_17 import table_4_8_dict_17_get
                # M_s = table_4_8_dict_17_get('URM', self.perf)
                # capacity_psi = (some base) / M_s   – need base value

                dc = demand_psi / capacity_psi
                self.results.append({
                    'check': f'Infill Panel Shear ({dirn}, storey {storey_idx})',
                    'demand_psi': demand_psi,
                    'capacity_psi': capacity_psi,
                    'dc_ratio': dc,
                    'pass': dc <= 1.0,
                    'note': 'Approximate capacity – verify with project data'
                })

    # ----------------------------------------------------------------------
    # (Optional) Helper for joint shear capacity – already defined earlier
    # ----------------------------------------------------------------------

    # ----------------------------------------------------------------------
    # Additional functions for type C2
    # ----------------------------------------------------------------------


    def _check_wall_redundancy(self):
        """Check that number of lines of shear walls >= 2 (Tables 17-24, 17-25)."""
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            
            lines_of_walls = d.get('lines_of_shear_walls', 0)
            required = 2
            passes = lines_of_walls >= required
            
            self.results.append({
                'check': f'Wall Redundancy ({dirn})',
                'key': 'REDUNDANCY',
                'tier2_ref': '5.5.1.1',
                'commentary_ref': 'A.3.2.1.1',
                'direction': dirn,
                'lines_of_walls': lines_of_walls,
                'required': required,
                'pass': passes,
                'note': '' if passes else f'Only {lines_of_walls} lines, need at least {required}'
            })


    def _check_wall_reinforcing(self):
        """
        Check wall reinforcing ratios and spacing.
        
        Table 17-24 (CP): vertical >= 0.0012, horizontal >= 0.0020
        Table 17-25 (IO): same ratios, plus spacing <= 18 in (457 mm)
        """
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            
            vertical_ratio = d.get('wall_vertical_reinf_ratio', 0)
            horizontal_ratio = d.get('wall_horizontal_reinf_ratio', 0)
            vertical_spacing = d.get('wall_vertical_bar_spacing_mm')
            horizontal_spacing = d.get('wall_horizontal_bar_spacing_mm')
            
            # Check ratios
            vertical_ok = vertical_ratio >= 0.0012
            horizontal_ok = horizontal_ratio >= 0.0020
            
            # Check spacing for IO
            spacing_ok = True
            issues = []
            
            if not vertical_ok:
                issues.append(f'Vertical ratio {vertical_ratio:.4f} < 0.0012')
            if not horizontal_ok:
                issues.append(f'Horizontal ratio {horizontal_ratio:.4f} < 0.0020')
            
            if self.perf.upper() == 'IO':
                max_spacing_mm = 457  # 18 inches
                if vertical_spacing and vertical_spacing > max_spacing_mm:
                    spacing_ok = False
                    issues.append(f'Vertical spacing {vertical_spacing}mm > 457mm')
                if horizontal_spacing and horizontal_spacing > max_spacing_mm:
                    spacing_ok = False
                    issues.append(f'Horizontal spacing {horizontal_spacing}mm > 457mm')
            
            passes = vertical_ok and horizontal_ok and spacing_ok
            
            self.results.append({
                'check': f'Wall Reinforcing ({dirn})',
                'key': 'REINFORCING_STEEL',
                'tier2_ref': '5.5.3.1.3',
                'commentary_ref': 'A.3.2.2.2',
                'direction': dirn,
                'vertical_ratio': vertical_ratio,
                'horizontal_ratio': horizontal_ratio,
                'vertical_spacing_mm': vertical_spacing,
                'horizontal_spacing_mm': horizontal_spacing,
                'pass': passes,
                'note': '; '.join(issues) if issues else 'Reinforcing ratios OK'
            })


    def _check_wall_aspect_ratios(self):
        """
        Check wall aspect ratio limits.
        
        Table 17-25 IO: aspect ratio < 4-to-1 for all walls.
        """
        walls = self.data.get('wall_information', [])
        
        if not walls:
            self.results.append({
                'check': 'Wall Aspect Ratio',
                'key': 'OVERTURNING',
                'tier2_ref': '5.5.3.1.4',
                'commentary_ref': 'A.3.2.2.4',
                'pass': None,
                'note': 'No wall information provided'
            })
            return
        
        for i, wall in enumerate(walls):
            # Calculate aspect ratio if not provided
            aspect = wall.get('aspect_ratio', 0)
            if aspect == 0 and wall.get('height_m') and wall.get('length_m'):
                aspect = wall['height_m'] / wall['length_m']
            
            if self.perf.upper() == 'IO':
                limit = 4.0
                passes = aspect < limit
                note = '' if passes else f'Aspect ratio {aspect:.2f} >= {limit}'
            else:
                # For CP, just report without pass/fail
                passes = None
                note = 'For information only (no CP limit)'
            
            self.results.append({
                'check': f'Wall Aspect Ratio (Wall {i+1}, {wall.get("orientation", "Unknown")})',
                'key': 'OVERTURNING',
                'tier2_ref': '5.5.3.1.4',
                'commentary_ref': 'A.3.2.2.4',
                'aspect_ratio': aspect,
                'pass': passes,
                'note': note
            })


    def _check_boundary_confinement(self):
        """
        Check boundary element confinement for walls with aspect ratio > 2.
        
        Table 17-25 IO: For walls with aspect ratio > 2-to-1,
        boundary elements confined with spirals or ties with spacing < 8db.
        """
        walls = self.data.get('wall_information', [])
        
        if not walls:
            return
        
        for i, wall in enumerate(walls):
            # Calculate aspect ratio if not provided
            aspect = wall.get('aspect_ratio', 0)
            if aspect == 0 and wall.get('height_m') and wall.get('length_m'):
                aspect = wall['height_m'] / wall['length_m']
            
            # Skip if aspect ratio <= 2
            if aspect <= 2.0:
                continue
            
            confined = wall.get('boundary_element_confined', False)
            tie_spacing = wall.get('confining_tie_spacing_mm')
            bar_size = wall.get('vertical_bar_size_mm')
            
            issues = []
            passes = True
            
            if not confined:
                passes = False
                issues.append('Boundary elements not confined')
            elif tie_spacing and bar_size:
                max_spacing = 8 * bar_size
                if tie_spacing > max_spacing:
                    passes = False
                    issues.append(f'Tie spacing {tie_spacing}mm > 8db ({max_spacing}mm)')
            
            self.results.append({
                'check': f'Boundary Confinement (Wall {i+1}, {wall.get("orientation", "Unknown")})',
                'key': 'CONFINEMENT_REINFORCING',
                'tier2_ref': '5.5.3.2.2',
                'commentary_ref': 'A.3.2.2.5',
                'aspect_ratio': aspect,
                'pass': passes if self.perf.upper() == 'IO' else None,
                'note': '; '.join(issues) if issues else 'Boundary confinement OK' if self.perf.upper() == 'IO' else 'Not required for CP'
            })


    def _check_wall_openings(self):
        """
        Check for trim reinforcement around wall openings.
        
        Table 17-25 IO: Added trim reinforcement required around all wall openings
        with dimension > 3 times wall thickness.
        """
        walls = self.data.get('wall_information', [])
        
        if not walls:
            return
        
        for i, wall in enumerate(walls):
            openings = wall.get('openings', [])
            if not openings:
                continue
            
            thickness = wall.get('thickness_m', 0.3)
            issues = []
            
            for j, opening in enumerate(openings):
                width = opening.get('width_m', 0)
                height = opening.get('height_m', 0)
                max_dim = max(width, height)
                
                if max_dim > 3 * thickness:
                    has_trim = opening.get('trim_reinforcing', False)
                    if not has_trim:
                        issues.append(f'Opening {j+1}: {max_dim:.2f}m > 3t ({3*thickness:.2f}m) missing trim')
            
            passes = len(issues) == 0
            
            self.results.append({
                'check': f'Wall Opening Reinforcing (Wall {i+1}, {wall.get("orientation", "Unknown")})',
                'key': 'WALL_REINFORCING_AT_OPENINGS',
                'tier2_ref': '5.5.3.1.5',
                'commentary_ref': 'A.3.2.2.6',
                'pass': passes if self.perf.upper() == 'IO' else None,
                'note': '; '.join(issues) if issues else 'Opening reinforcing OK' if passes else ''
            })


    def _check_wall_thickness(self):
        """
        Check minimum wall thickness requirements.
        
        Table 17-25 IO: Thickness not less than 1/25 the unsupported height or length,
        nor less than 4 in (101 mm).
        """
        walls = self.data.get('wall_information', [])
        storey_heights = self.data.get('storey_elevations_m', [])
        
        if not walls:
            return
        
        for i, wall in enumerate(walls):
            storey = wall.get('storey', 0)
            
            # Determine unsupported height
            if storey < len(storey_heights):
                if storey == 0:
                    # First storey - height from ground to first floor
                    height = storey_heights[0]
                else:
                    # Upper storeys
                    height = storey_heights[storey-1] - storey_heights[storey]
            else:
                height = wall.get('height_m', 4.0)
            
            thickness = wall.get('thickness_m', 0.3)
            length = wall.get('length_m', 0)
            
            # Calculate minimum required thickness
            min_by_height = height / 25
            min_by_length = length / 25 if length > 0 else 0
            required_min_m = max(min_by_height, min_by_length, 0.101)  # 101 mm = 0.101 m
            
            passes = thickness >= required_min_m
            
            self.results.append({
                'check': f'Wall Thickness (Wall {i+1}, {wall.get("orientation", "Unknown")})',
                'key': 'WALL_THICKNESS',
                'tier2_ref': '5.5.3.1.2',
                'commentary_ref': 'A.3.2.2.7',
                'thickness_mm': thickness * 1000,
                'required_min_mm': required_min_m * 1000,
                'pass': passes if self.perf.upper() == 'IO' else None,
                'note': '' if passes else f'Thickness {thickness*1000:.0f}mm < required {required_min_m*1000:.0f}mm'
            })


    def _check_coupling_beams(self):
        """
        Check coupling beam requirements.
        
        Table 17-25 IO: Coupling beams have capacity to develop uplift capacity of adjacent wall.
        Table 17-24 CP: Ends of both walls supported to resist vertical loads.
        """
        for dirn in ('X', 'Y'):
            d = self.data['directions'].get(dirn)
            if not d:
                continue
            
            coupling_present = d.get('coupling_beams_present', False)
            
            if not coupling_present:
                continue
            
            details = d.get('coupling_beam_details', {})
            issues = []
            
            # Check if walls are supported (CP requirement)
            walls_supported = details.get('walls_supported', True)
            if not walls_supported and self.perf.upper() == 'CP':
                issues.append('Walls not supported at ends')
            
            # Check capacity (IO requirement)
            if self.perf.upper() == 'IO':
                span = details.get('span_m', 0)
                depth = details.get('depth_m', 0)
                reinf_ratio = details.get('reinforcement_ratio', 0)
                
                if depth < span / 4:
                    issues.append(f'Beam depth {depth:.2f}m may be inadequate for span {span:.2f}m')
                if reinf_ratio < 0.01:
                    issues.append(f'Reinforcement ratio {reinf_ratio:.4f} may be low')
            
            passes = len(issues) == 0
            
            self.results.append({
                'check': f'Coupling Beams ({dirn})',
                'key': 'COUPLING_BEAMS',
                'tier2_ref': '5.5.3.2.1',
                'commentary_ref': 'A.3.2.2.3',
                'direction': dirn,
                'pass': passes,
                'note': '; '.join(issues) if issues else 'Coupling beams OK'
            })


    def _check_secondary_frames(self):
        """
        Check that secondary frames form a complete vertical load system.
        """
        frame_system = self.data.get('frame_system', {})
        present = frame_system.get('secondary_frames_present', False)
        
        if not present:
            self.results.append({
                'check': 'Secondary Frames',
                'key': 'COMPLETE_FRAMES',
                'tier2_ref': '5.5.2.5.1',
                'commentary_ref': 'A.3.1.6.1',
                'pass': True,
                'note': 'No secondary frames present'
            })
            return
        
        complete_system = frame_system.get('complete_vertical_system', True)
        
        self.results.append({
            'check': 'Secondary Frames',
            'key': 'COMPLETE_FRAMES',
            'tier2_ref': '5.5.2.5.1',
            'commentary_ref': 'A.3.1.6.1',
            'pass': complete_system,
            'note': '' if complete_system else 'Secondary frames do not form complete vertical load system'
        })


    def _check_diaphragm_continuity(self):
        """
        Check diaphragm continuity requirements.
        """
        dia = self.data.get('diaphragm_info', {})
        
        continuous = dia.get('continuous', True)
        split_level = dia.get('split_level', False)
        expansion_joints = dia.get('expansion_joints', False)
        
        passes = continuous and not split_level and not expansion_joints
        
        issues = []
        if not continuous:
            issues.append('Diaphragm not continuous')
        if split_level:
            issues.append('Split-level floors present')
        if expansion_joints:
            issues.append('Expansion joints present')
        
        self.results.append({
            'check': 'Diaphragm Continuity',
            'key': 'DIAPHRAGM_CONTINUITY',
            'tier2_ref': '5.6.1.1',
            'commentary_ref': 'A.4.1.1',
            'pass': passes,
            'note': '; '.join(issues) if issues else 'Diaphragm continuous'
        })


    def _check_diaphragm_openings_at_walls(self):
        """
        Check diaphragm openings adjacent to shear walls.
        
        Table 17-24 CP: Openings < 25% of wall length
        Table 17-25 IO: Openings < 15% of wall length
        """
        dia = self.data.get('diaphragm_info', {})
        openings = dia.get('openings', [])
        
        if not openings:
            self.results.append({
                'check': 'Diaphragm Openings at Walls',
                'key': 'OPENINGS_AT_SHEAR_WALLS',
                'tier2_ref': '5.6.1.3',
                'commentary_ref': 'A.4.1.4',
                'pass': True,
                'note': 'No diaphragm openings'
            })
            return
        
        limit = 0.15 if self.perf.upper() == 'IO' else 0.25
        issues = []
        
        for i, opening in enumerate(openings):
            if opening.get('adjacent_to_wall', False):
                size_percent = opening.get('size_percent', 0) / 100
                wall_length = opening.get('wall_length_affected_m', 0)
                
                if size_percent > limit:
                    issues.append(f'Opening {i+1}: {size_percent*100:.0f}% > {limit*100:.0f}% limit')
        
        passes = len(issues) == 0
        
        self.results.append({
            'check': 'Diaphragm Openings at Walls',
            'key': 'OPENINGS_AT_SHEAR_WALLS',
            'tier2_ref': '5.6.1.3',
            'commentary_ref': 'A.4.1.4',
            'pass': passes,
            'note': '; '.join(issues) if issues else 'Opening sizes within limits'
        })


    def _check_diaphragm_cross_ties(self):
        """
        Check for continuous cross ties between diaphragm chords.
        """
        dia = self.data.get('diaphragm_info', {})
        present = dia.get('cross_ties_present', False)
        
        self.results.append({
            'check': 'Diaphragm Cross Ties',
            'key': 'CROSS_TIES',
            'tier2_ref': '5.6.1.2',
            'commentary_ref': 'A.4.1.2',
            'pass': present,
            'note': '' if present else 'Continuous cross ties not present'
        })


    def _check_diaphragm_sheathing(self):
        """
        Check diaphragm sheathing requirements.
        """
        dia = self.data.get('diaphragm_info', {})
        aspect = dia.get('aspect_ratio', 1.0)
        span_m = dia.get('span_m', 0)
        span_ft = m2ft(span_m)
        
        # Check straight sheathing
        if dia.get('straight_sheathing', False):
            limit = 1.0 if self.perf.upper() == 'IO' else 2.0
            passes = aspect < limit
            
            self.results.append({
                'check': 'Straight Sheathing',
                'key': 'STRAIGHT_SHEATHING',
                'tier2_ref': '5.6.2',
                'commentary_ref': 'A.4.2.1',
                'aspect_ratio': aspect,
                'limit': limit,
                'pass': passes,
                'note': '' if passes else f'Aspect ratio {aspect:.2f} >= {limit}'
            })
        
        # Check diagonal/unblocked sheathing
        if dia.get('diagonal_sheathing', False) or dia.get('unblocked', False):
            span_limit_ft = 30.0 if self.perf.upper() == 'IO' else 40.0
            aspect_limit = 3.0 if self.perf.upper() == 'IO' else 4.0
            
            span_ok = span_ft < span_limit_ft
            aspect_ok = aspect <= aspect_limit
            
            passes = span_ok and aspect_ok
            issues = []
            if not span_ok:
                issues.append(f'Span {span_ft:.1f}ft > {span_limit_ft}ft')
            if not aspect_ok:
                issues.append(f'Aspect ratio {aspect:.2f} > {aspect_limit}')
            
            self.results.append({
                'check': 'Diagonal/Unblocked Sheathing',
                'key': 'DIAGONALLY_SHEATHED',
                'tier2_ref': '5.6.2',
                'commentary_ref': 'A.4.2.3',
                'span_ft': span_ft,
                'aspect_ratio': aspect,
                'pass': passes,
                'note': '; '.join(issues) if issues else 'Sheathing OK'
            })


    def _check_metal_deck(self):
        """
        Check untopped metal deck diaphragms.
        
        Table 17-25: Spans < 40 ft (12.2 m), aspect ratio < 4-to-1.
        """
        dia = self.data.get('diaphragm_info', {})
        
        if not dia.get('untopped_metal_deck', False):
            return
        
        span_m = dia.get('span_m', 0)
        span_ft = m2ft(span_m)
        aspect = dia.get('aspect_ratio', 1.0)
        
        span_ok = span_ft < 40.0
        aspect_ok = aspect < 4.0
        
        passes = span_ok and aspect_ok
        issues = []
        if not span_ok:
            issues.append(f'Span {span_ft:.1f}ft >= 40ft')
        if not aspect_ok:
            issues.append(f'Aspect ratio {aspect:.2f} >= 4.0')
        
        self.results.append({
            'check': 'Untopped Metal Deck',
            'key': 'NONCONCRETE_FILLED',
            'tier2_ref': '5.6.3',
            'commentary_ref': 'A.4.3.1',
            'span_ft': span_ft,
            'aspect_ratio': aspect,
            'pass': passes,
            'note': '; '.join(issues) if issues else 'Metal deck OK'
        })


    def _check_foundation_dowels(self):
        """
        Check foundation dowels for walls.
        """
        found = self.data.get('foundation_info', {})
        dowels = found.get('wall_dowels', {})
        
        if not dowels:
            self.results.append({
                'check': 'Foundation Dowels',
                'key': 'FOUNDATION_DOWELS',
                'tier2_ref': '5.7.3.4',
                'commentary_ref': 'A.5.3.5',
                'pass': None,
                'note': 'No dowel information provided'
            })
            return
        
        size_equal = dowels.get('size_equal_to_wall_reinf', False)
        spacing_equal = dowels.get('spacing_equal_to_wall_reinf', False)
        develop_strength = dowels.get('able_to_develop_strength', False)
        
        passes = size_equal and spacing_equal and develop_strength
        issues = []
        if not size_equal:
            issues.append('Dowel size not equal to wall reinforcement')
        if not spacing_equal:
            issues.append('Dowel spacing not equal to wall reinforcement')
        if not develop_strength:
            issues.append('Dowels unable to develop required strength')
        
        self.results.append({
            'check': 'Foundation Dowels',
            'key': 'FOUNDATION_DOWELS',
            'tier2_ref': '5.7.3.4',
            'commentary_ref': 'A.5.3.5',
            'pass': passes,
            'note': '; '.join(issues) if issues else 'Dowels OK'
        })


    def _check_pile_cap_uplift(self):
        """
        Check pile cap uplift resistance.
        
        Table 17-25 IO: Pile caps have top reinforcement, piles anchored,
        and can develop tensile capacity.
        Table 17-24 CP: Only require top reinforcement and anchored piles.
        """
        found = self.data.get('foundation_info', {})
        pile_caps = found.get('pile_caps', {})
        
        if not pile_caps:
            return
        
        top_reinf = pile_caps.get('top_reinforcement', False)
        piles_anchored = pile_caps.get('piles_anchored', False)
        tensile_capacity = pile_caps.get('tensile_capacity_developed', False)
        
        passes = top_reinf and piles_anchored
        issues = []
        
        if not top_reinf:
            issues.append('No top reinforcement')
        if not piles_anchored:
            issues.append('Piles not anchored')
        
        if self.perf.upper() == 'IO':
            passes = passes and tensile_capacity
            if not tensile_capacity:
                issues.append('Cannot develop tensile capacity')
        
        self.results.append({
            'check': 'Pile Cap Uplift',
            'key': 'UPLIFT_AT_PILE_CAPS',
            'tier2_ref': '5.7.3.5',
            'commentary_ref': 'A.5.3.8',
            'pass': passes,
            'note': '; '.join(issues) if issues else 'Pile caps OK'
        })




# ============================================================================
# ASCE 41-17 Tier 1 Checklists
# ============================================================================

# Table 17-2: Collapse Prevention Basic Configuration Checklist
# (For Low, Moderate, and High Seismicity)
CP_BASIC_CONFIGURATION_CHECKLIST = {
    # Very low seismicity
    'LOAD_PATH': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The structure contains a complete, well-defined load path, including structural elements and connections, that serves to transfer the inertial forces associated with the mass of all elements of the building to the foundation.',
        'tier2_ref': '5.4.1.1',
        'commentary_ref': 'A.2.1.1'
    },
    'WALL_ANCHORAGE': {
        'seismicity': ['Very Low'],
        'statement': 'Exterior concrete or masonry walls that are dependent on the diaphragm for lateral support are anchored for out-of-plane forces at each diaphragm level with steel anchors, reinforcing dowels, or straps that are developed into the diaphragm. Connections have adequate strength to resist the connection force calculated in the Quick Check procedure of Section 4.4.3.7.',
        'tier2_ref': '5.7.1.1',
        'commentary_ref': 'A.5.1.1'
    },
    # Low Seismicity and above
    'ADJACENT_BUILDINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The clear distance between the building being evaluated and any adjacent building is greater than 0.25% of the height of the shorter building in low seismicity, 0.5% in moderate seismicity, and 1.5% in high seismicity.',
        'tier2_ref': '5.4.1.2',
        'commentary_ref': 'A.2.1.2'
    },
    'MEZZANINES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Interior mezzanine levels are braced independently from the main structure or are anchored to the seismic-force-resisting elements of the main structure.',
        'tier2_ref': '5.4.1.3',
        'commentary_ref': 'A.2.1.3'
    },
    'WEAK_STORY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The sum of the shear strengths of the seismic-force-resisting system in any story in each direction is not less than 80% of the strength in the adjacent story above.',
        'tier2_ref': '5.4.2.1',
        'commentary_ref': 'A.2.2.2'
    },
    'SOFT_STORY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The stiffness of the seismic-force-resisting system in any story is not less than 70% of the seismic-force-resisting system stiffness in an adjacent story above or less than 80% of the average seismic-force-resisting system stiffness of the three stories above.',
        'tier2_ref': '5.4.2.2',
        'commentary_ref': 'A.2.2.3'
    },
    'VERTICAL_IRREGULARITIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All vertical elements in the seismic-force-resisting system are continuous to the foundation.',
        'tier2_ref': '5.4.2.3',
        'commentary_ref': 'A.2.2.4'
    },
    'GEOMETRY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There are no changes in the net horizontal dimension of the seismic-force-resisting system of more than 30% in a story relative to adjacent stories, excluding one-story penthouses and mezzanines.',
        'tier2_ref': '5.4.2.4',
        'commentary_ref': 'A.2.2.5'
    },
    'MASS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is no change in effective mass of more than 50% from one story to the next. Light roofs, penthouses, and mezzanines need not be considered.',
        'tier2_ref': '5.4.2.5',
        'commentary_ref': 'A.2.2.6'
    },
    'TORSION': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The estimated distance between the story center of mass and the story center of rigidity is less than 20% of the building width in either plan dimension.',
        'tier2_ref': '5.4.2.6',
        'commentary_ref': 'A.2.2.7'
    },
    # ----- Low Seismicity (Additional Items) -----
    'LIQUEFACTION': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'Liquefaction-susceptible, saturated, loose granular soils that could jeopardize the building\'s seismic performance do not exist in the foundation soils at depths within 50 ft (15.2 m) under the building.',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.1'
    },
    'SLOPE_FAILURE': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'The building site is located away from potential earthquake-induced slope failures or rockfalls so that it is unaffected by such failures or is capable of accommodating any predicted movements without failure.',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.2'
    },
    'SURFACE_FAULT_RUPTURE': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'Surface fault rupture and surface displacement at the building site are not anticipated.',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.3'
    },
    
    # ----- Moderate and High Seismicity (Additional Items) -----
    'OVERTURNING': {
        'seismicity': ['High'],
        'statement': 'The ratio of the least horizontal dimension of the seismic-force-resisting system at the foundation level to the building height (base-height) is greater than 0.6Sa.',
        'tier2_ref': '5.4.3.3',
        'commentary_ref': 'A.6.2.1'
    },
    'TIES_BETWEEN_FOUNDATION_ELEMENTS': {
        'seismicity': ['High'],
        'statement': 'The foundation has ties adequate to resist seismic forces where footings, piles, and piers are not restrained by beams, slabs, or soils classified as Site Class A, B, or C.',
        'tier2_ref': '5.4.3.4',
        'commentary_ref': 'A.6.2.2'
    }
}

IO_BASIC_CONFIGURATION_CHECKLIST = {
    # ----- Building System—General (Very Low Seismicity) -----
    # ----- Building System—General (Low Seismicity and above) -----
    'LOAD_PATH': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The structure contains a complete, well-defined load path, including structural elements and connections, that serves to transfer the inertial forces associated with the mass of all elements of the building to the foundation.',
        'tier2_ref': '5.4.1.1',
        'commentary_ref': 'A.2.1.1'
    },
    'ADJACENT_BUILDINGS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The clear distance between the building being evaluated and any adjacent building is greater than 0.5% of the height of the shorter building in low seismicity, 1.0% in moderate seismicity, and 2.0% in high seismicity.',
        'tier2_ref': '5.4.1.2',
        'commentary_ref': 'A.2.1.2'
    },
    'MEZZANINES': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Interior mezzanine levels are braced independently from the main structure and are anchored to the seismic-force-resisting elements of the main structure with connections capable of developing the full capacity of the mezzanine framing.',
        'tier2_ref': '5.4.1.3',
        'commentary_ref': 'A.2.1.3'
    },
    
    # ----- Building System—Building Configuration (Low Seismicity and above) -----
    'WEAK_STORY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The sum of the shear strengths of the seismic-force-resisting system in any story in each direction is not less than 90% of the strength in the adjacent story above.',
        'tier2_ref': '5.4.2.1',
        'commentary_ref': 'A.2.2.2'
    },
    'SOFT_STORY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The stiffness of the seismic-force-resisting system in any story is not less than 80% of the seismic-force-resisting system stiffness in an adjacent story above or less than 90% of the average seismic-force-resisting system stiffness of the three stories above.',
        'tier2_ref': '5.4.2.2',
        'commentary_ref': 'A.2.2.3'
    },
    'VERTICAL_IRREGULARITIES': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'All vertical elements in the seismic-force-resisting system are continuous to the foundation and have positive connection at all splices and connections.',
        'tier2_ref': '5.4.2.3',
        'commentary_ref': 'A.2.2.4'
    },
    'GEOMETRY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'There are no changes in the net horizontal dimension of the seismic-force-resisting system of more than 20% in a story relative to adjacent stories, excluding one-story penthouses and mezzanines.',
        'tier2_ref': '5.4.2.4',
        'commentary_ref': 'A.2.2.5'
    },
    'MASS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'There is no change in effective mass of more than 40% from one story to the next. Light roofs, penthouses, and mezzanines need not be considered.',
        'tier2_ref': '5.4.2.5',
        'commentary_ref': 'A.2.2.6'
    },
    'TORSION': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The estimated distance between the story center of mass and the story center of rigidity is less than 15% of the building width in either plan dimension.',
        'tier2_ref': '5.4.2.6',
        'commentary_ref': 'A.2.2.7'
    },
    
    # ----- Low Seismicity (Additional Items) -----
    'LIQUEFACTION': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Liquefaction-susceptible, saturated, loose granular soils that could jeopardize the building\'s seismic performance do not exist in the foundation soils at depths within 75 ft (22.9 m) under the building. If present, the soils are stabilized or the foundation is designed to accommodate liquefaction effects.',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.1'
    },
    'SLOPE_FAILURE': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The building site is located away from potential earthquake-induced slope failures or rockfalls, or the building is designed to accommodate the effects of such failures without damage.',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.2'
    },
    'SURFACE_FAULT_RUPTURE': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Surface fault rupture and surface displacement at the building site are not anticipated (? or the building is designed to accommodate the effects of such displacement).',
        'tier2_ref': '5.4.3.1',
        'commentary_ref': 'A.6.1.3'
    },
    
    # ----- Moderate and High Seismicity (Additional Items) -----
    'OVERTURNING': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'The ratio of the least horizontal dimension of the seismic-force-resisting system at the foundation level to the building height (base-height) is greater than 0.8Sa.',
        'tier2_ref': '5.4.3.3',
        'commentary_ref': 'A.6.2.1'
    },
    'TIES_BETWEEN_FOUNDATION_ELEMENTS': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'The foundation has ties capable of developing the full seismic force demands where footings, piles, and piers are not restrained by beams, slabs, or soils classified as Site Class A, B, or C.',
        'tier2_ref': '5.4.3.4',
        'commentary_ref': 'A.6.2.2'
    },

}

# ============================================================================
# W1: WOOD LIGHT FRAMES
# ============================================================================

# ...

# ============================================================================
# C1: CONCRETE MOMENT FRAMES
# ============================================================================

# ============================================================================
# ASCE 41-17 Tier 1 Checklists
# ============================================================================

# ============================================================================
# C1: CONCRETE MOMENT FRAMES
# ============================================================================

# Table 17-22: Collapse Prevention Structural Checklist for Building Type C1
C1_CP_CHECKLIST = {
    # Low Seismicity
    'REDUNDANCY_1': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The number of lines of moment frames in each principal direction is greater than or equal to 2. The number of bays of moment frames in each line is greater than or equal to 3.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.1.1.1'
    },
    'COLUMN_AXIAL_STRESS_CHECK': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The axial stress caused by unfactored gravity loads in columns subjected to overturning forces because of seismic demands is less than 0.20f'c. Alternatively, the axial stress caused by overturning forces alone, calculated using the Quick Check procedure of Section 4.4.3.6, is less than 0.30f\'c.",
        'tier2_ref': '5.5.2.1.3',
        'commentary_ref': 'A.3.1.4.2'
    },
    'CONCRETE_COLUMNS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "All concrete columns are doweled into the foundation with a minimum of four bars",
        'tier2_ref': '5.7.3.1',
        'commentary_ref': 'A.5.3.2'
    },
    # Moderate Seismicity
    'REDUNDANCY_2': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'The number of bays of moment frames in each line is greater than or equal to 2. The number of bays of moment frames in each line is greater than or equal to 3.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.1.1.1'
    },
    'INTERFERING_WALLS': {
        'seismicity': ['Moderate', 'High'],
        'statement': 'All concrete and masonry infill walls placed in moment frames are isolated from the structural members.',
        'tier2_ref': '5.5.2.1.1',
        'commentary_ref': 'A.3.1.2.1'
    },
    'COLUMN_SHEAR_STRESS_CHECK': {
        'seismicity': ['Moderate', 'High'],
        'statement': "The shear stress in the concrete moment frame columns, calculated using the Quick Check procedure of Section 4.4.3.2, is less than the greater of 100 lb/in² (0.69 MPa) or 2√f\'c.",
        'tier2_ref': '5.5.2.1.4',
        'commentary_ref': 'A.3.1.4.1'
    },
    'FLAT_SLAB_FRAMES': {
        'seismicity': ['Moderate', 'High'],
        'statement': "The seismic-force-resisting system is not a frame consisting of columns and a flat slab or plate without beams.",
        'tier2_ref': '5.5.2.3.1',
        'commentary_ref': 'A.3.1.4.3'
    },
    # High Seismicity (additional items)
    'PRESTRESSED_FRAME_ELEMENTS': {
        'seismicity': ['High'],
        'statement': "The seismic-force-resisting frames do not include any prestressed or post-tensioned elements where the average prestress exceeds the lesser of 700 lb/in.2 (4.83 MPa) or f’c /6 at potential hinge locations. The average prestress is calculated in accordance with the Quick Check procedure of Section 4.4.3.8.",
        'tier2_ref': '5.5.2.3.2',
        'commentary_ref': 'A.3.1.4.4'
    },
    'CAPTIVE_COLUMNS':{
        'seismicity': ['High'],
        'statement': "There are no columns at a level with height/depth ratios less than 50% of the nominal height/depth ratio of the typical columns at that level.",
        'tier2_ref': '5.5.2.3.3',
        'commentary_ref': 'A.3.1.4.5'
    },
    'NO_SHEAR_FAILURES': {
        'seismicity': ['High'],
        'statement': "The shear capacity of frame members is able to develop the moment capacity at the ends of the members.",
        'tier2_ref': '5.5.2.3.4',
        'commentary_ref': 'A.3.1.4.6'
    },
    'STRONG_COLUMN_WEAK_BEAM': {
        'seismicity': ['High'],
        'statement': 'The sum of the moment capacity of the columns is 20% greater than that of the beams at frame joints.',
        'tier2_ref': '5.5.2.1.5',
        'commentary_ref': 'A.3.1.4.7'
    },
    'BEAM_BARS':{
        'seismicity': ['High'],
        'statement': "At least two longitudinal top and two longitudinal bottom bars extend continuously throughout the length of each frame beam. At least 25% of the longitudinal bars provided at the joints for either positive or negative moment are continuous throughout the length of the members.",
        'tier2_ref': '5.5.2.3.5',
        'commentary_ref': 'A.3.1.4.8'
        },
    'COLUMN_BAR_SPLICES': {
        'seismicity': ['High'],
        'statement': 'All column-bar lap splice lengths are greater than 35db and are enclosed by ties spaced at or less than 8db. Alternatively, column bars are spliced with mechanical couplers with a capacity of at least 1.25 times the nominal yield strength of the spliced bar.',
        'tier2_ref': '5.5.2.3.6',
        'commentary_ref': 'A.3.1.4.9'
    },
    'BEAM_BAR_SPLICES': {
        'seismicity': ['High'],
        'statement': "The lap splices or mechanical couplers for longitudinal beam reinforcing are not located within lb /4 of the joints and are not located in the vicinity of potential plastic hinge locations.",
        'tier2_ref': '5.5.2.3.6',
        'commentary_ref': 'A.3.1.4.10'
    },
    'COLUMN_TIE_SPACING': {
        'seismicity': ['High'],
        'statement': "Frame columns have ties spaced at or less than d/4 throughout their length and at or less than 8db at all potential plastic hinge locations.",
        'tier2_ref': '5.5.2.3.7',
        'commentary_ref': 'A.3.1.4.11'
    },
    'STIRRUP_SPACING': {
        'seismicity': ['High'],
        'statement': 'All beams have stirrups spaced at or less than d/2 throughout their length. At potential plastic hinge locations, stirrups are spaced at or less than the minimum of 8db or d/4.',
        'tier2_ref': '5.5.2.3.7',
        'commentary_ref': 'A.3.1.4.12'
    },
    'JOINT_TRANSVERSE_REINFORCING': {
        'seismicity': ['High'],
        'statement': 'Beam–column joints have ties spaced at or less than 8db.',
        'tier2_ref': '5.5.2.3.8',
        'commentary_ref': 'A.3.1.4.13'
    },
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['High'],
        'statement': 'Flat slabs or plates not part of the seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'DIAPHRAGM_CONTINUITY':{
        'seismicity': ['High'],
        'statement': "The diaphragms are not composed of split-level floors and do not have expansion joints.",
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'UPLIFT_AT_PILE_CAPS':{
        'seismicity': ['High'],
        'statement': "Pile caps have top reinforcement, and piles are anchored to the pile caps.",
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    },
}

# 'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],

# Table 17-23: Immediate Occupancy Structural Checklist for Building Type C1
C1_IO_CHECKLIST = {
    # Very Low Seismicity (and also apply to Low/Moderate unless superseded)
    # Seismic Force-resisting system
    'REDUNDANCY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The number of lines of moment frames in each principal direction is greater than or equal to 2. The number of bays of moment frames in each line is greater than or equal to 3.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.1.1.1'
    },
    'INTERFERING_WALLS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'All concrete and masonry infill walls placed in moment frames are isolated from structural members.',
        'tier2_ref': '5.5.2.1.1',
        'commentary_ref': 'A.3.1.2.1'
    },
    'COLUMN_SHEAR_STRESS_CHECK': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': "The shear stress in the concrete moment frame columns, calculated using the Quick Check procedure of Section 4.4.3.2, is less than the greater of 100 lb/in² (0.69 MPa) or 2√f\'c.",
        'tier2_ref': '5.5.2.1.4',
        'commentary_ref': 'A.3.1.4.1'
    },
    'COLUMN_AXIAL_STRESS_CHECK': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': "The axial stress caused by unfactored gravity loads in columns subjected to overturning demands is less than 0.13f'c. Alternatively, the axial stress caused by overturning forces alone, calculated using the Quick Check procedure of Section 4.4.3.6, is less than 0.30f'c.",
        'tier2_ref': '5.5.2.1.3',
        'commentary_ref': 'A.3.1.4.2'
    },
    # Connections
    'CONCRETE_COLUMNS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': "All concrete columns are doweled into the foundation, and the dowels are able to develop the tensile capacity of reinforcement in columns of the seismic-force-resisting system.",
        'tier2_ref': '5.7.3.1',
        'commentary_ref': 'A.5.3.2'
    },
    
    # Low, Moderate, and High Seismicity (more restrictive items)
	'FLAT_SLAB_FRAMES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The seismic-force-resisting system is not a frame consisting of columns and a flat slab or plate without beams.",
        'tier2_ref': '5.5.2.3.1',
        'commentary_ref': 'A.3.1.4.3'
    },
    'PRESTRESSED_FRAME_ELEMENTS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The seismic-force-resisting frames do not include any prestressed or post-tensioned elements where the average prestress exceeds the lesser of 700 lb/in.2 (4.83 MPa) or f’c /6 at potential hinge locations. The average prestress is calculated in accordance with the Quick Check procedure of Section 4.4.3.8.",
        'tier2_ref': '5.5.2.3.2',
        'commentary_ref': 'A.3.1.4.4'
    },
    'CAPTIVE_COLUMNS':{
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "There are no columns at a level with height/depth ratios less than 50% of the nominal height/depth ratio of the typical columns at that level.",
        'tier2_ref': '5.5.2.3.3',
        'commentary_ref': 'A.3.1.4.5'
    },
    'NO_SHEAR_FAILURES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The shear capacity of frame members is able to develop the moment capacity at the ends of the members.",
        'tier2_ref': '5.5.2.3.4',
        'commentary_ref': 'A.3.1.4.6'
    },
    'STRONG_COLUMN_WEAK_BEAM': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The sum of the moment capacity of the columns is 20% greater than that of the beams at frame joints.',
        'tier2_ref': '5.5.2.1.5',
        'commentary_ref': 'A.3.1.4.7'
    },
    'BEAM_BARS':{
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "At least two longitudinal top and two longitudinal bottom bars extend continuously throughout the length of each frame beam. At least 25% of the longitudinal bars provided at the joints for either positive or negative moment are continuous throughout the length of the members.",
        'tier2_ref': '5.5.2.3.5',
        'commentary_ref': 'A.3.1.4.8'
        },
    'COLUMN_BAR_SPLICES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All column-bar lap splice lengths are greater than 35db and are enclosed by ties spaced at or less than 8db. Alternatively, column bars are spliced with mechanical couplers with a capacity of at least 1.25 times the nominal yield strength of the spliced bar.',
        'tier2_ref': '5.5.2.3.6',
        'commentary_ref': 'A.3.1.4.9'
    },
    'BEAM_BAR_SPLICES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The lap splices or mechanical couplers for longitudinal beam reinforcing are not located within lb /4 of the joints and are not located in the vicinity of potential plastic hinge locations.",
        'tier2_ref': '5.5.2.3.6',
        'commentary_ref': 'A.3.1.4.10'
    },
    'COLUMN_TIE_SPACING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "Frame columns have ties spaced at or less than d/4 throughout their length and at or less than 8db at all potential plastic hinge locations.",
        'tier2_ref': '5.5.2.3.7',
        'commentary_ref': 'A.3.1.4.11'
    },
    'STIRRUP_SPACING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All beams have stirrups spaced at or less than d/2 throughout their length. At potential plastic hinge locations, stirrups are spaced at or less than the minimum of 8db or d/4.',
        'tier2_ref': '5.5.2.3.7',
        'commentary_ref': 'A.3.1.4.12'
    },
    'JOINT_TRANSVERSE_REINFORCING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Beam–column joints have ties spaced at or less than 8db.',
        'tier2_ref': '5.5.2.3.8',
        'commentary_ref': 'A.3.1.4.13'
    },
    'JOINT_ECCENTRICITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There are no eccentricities larger than 20% of the smallest column plan dimension between girder and column centerlines.',
        'tier2_ref': '5.5.2.3.9',
        'commentary_ref': 'A.3.1.4.14'
    },
    'STIRRUP_AND_TIE_HOOKS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The beam stirrups and column ties are anchored into the member cores with hooks of 135 degrees or more.',
        'tier2_ref': '5.5.2.3.10',
        'commentary_ref': 'A.3.1.4.15'
    },
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components and are Compliant with the following items: COLUMN-BAR SPLICES, BEAM-BAR SPLICES, COLUMN-TIE SPACING, STIRRUP SPACING, and STIRRUP AND TIE HOOKS.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Flat slabs or plates not part of the seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'DIAPHRAGM_CONTINUITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The diaphragms are not composed of split-level floors and do not have expansion joints.',
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'PLAN_IRREGULARITIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is tensile capacity to develop the strength of the diaphragm at reentrant corners or other locations of plan irregularities.',
        'tier2_ref': '5.6.1.4',
        'commentary_ref': 'A.4.1.7'
    },
    'DIAPHRAGM_REINFORCEMENT_AT_OPENINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is reinforcing around all diaphragm openings larger than 50% of the building width in either major plan dimension.',
        'tier2_ref': '5.6.1.5',
        'commentary_ref': 'A.4.1.8'
    },
    'UPLIFT_AT_PILE_CAPS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Pile caps have top reinforcement, and piles are anchored to the pile caps; the pile cap reinforcement and pile anchorage are able to develop the tensile capacity of the piles.',
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    }
}


# ============================================================================
# C2: CONCRETE SHEAR WALLS WITH STIFF DIAPHRAGMS
# C2A: CONCRETE SHEAR WALLS WITH FLEXIBLE DIAPHRAGMS
# ============================================================================

# Table 17-24: Collapse Prevention Structural Checklist for Building Types C2 and C2a
C2_CP_CHECKLIST = {
    # Low and Moderate Seismicity
    'COMPLETE_FRAMES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Steel or concrete frames classified as secondary components form a complete vertical-load-carrying system.',
        'tier2_ref': '5.5.2.5.1',
        'commentary_ref': 'A.3.1.6.1'
    },
    'REDUNDANCY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The number of lines of shear walls in each principal direction is greater than or equal to 2.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.2.1.1'
    },
    'SHEAR_STRESS_CHECK': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The shear stress in the concrete shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than the greater of 100 lb/in² (0.69 MPa) or 2√f\'c.",
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.2.1'
    },
    'REINFORCING_STEEL': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The ratio of reinforcing steel area to gross concrete area is not less than 0.0012 in the vertical direction and 0.0020 in the horizontal direction.',
        'tier2_ref': '5.5.3.1.3',
        'commentary_ref': 'A.3.2.2.2'
    },
    'WALL_ANCHORAGE_FLEXIBLE': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Exterior concrete or masonry walls that are dependent on flexible diaphragms for lateral support are anchored for out-of-plane forces at each diaphragm level with steel anchors, reinforcing dowels, or straps that are developed into the diaphragm. Connections have strength to resist the connection force calculated in the Quick Check procedure of Section 4.4.3.7.',
        'tier2_ref': '5.7.1.1',
        'commentary_ref': 'A.5.1.1',
        'diaphragm_type': 'Flexible'
    },
    'TRANSFER_TO_SHEAR_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragms are connected for transfer of seismic forces to the shear walls.',
        'tier2_ref': '5.7.2',
        'commentary_ref': 'A.5.2.1'
    },
    'FOUNDATION_DOWELS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Wall reinforcement is doweled into the foundation with vertical bars equal in size and spacing to the vertical wall reinforcing directly above the foundation.',
        'tier2_ref': '5.7.3.4',
        'commentary_ref': 'A.5.3.5'
    },
    
    # High Seismicity (additional items)
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['High'],
        'statement': 'Flat slabs or plates not part of the seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'COUPLING_BEAMS': {
        'seismicity': ['High'],
        'statement': 'The ends of both walls to which the coupling beam is attached are supported at each end to resist vertical loads caused by overturning.',
        'tier2_ref': '5.5.3.2.1',
        'commentary_ref': 'A.3.2.2.3'
    },
    'DIAPHRAGM_CONTINUITY': {
        'seismicity': ['High'],
        'statement': 'The diaphragms are not composed of split-level floors and do not have expansion joints.',
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'OPENINGS_AT_SHEAR_WALLS': {
        'seismicity': ['High'],
        'statement': 'Diaphragm openings immediately adjacent to the shear walls are less than 25% of the wall length.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.4'
    },
    'CROSS_TIES': {
        'seismicity': ['High'],
        'statement': 'There are continuous cross ties between diaphragm chords.',
        'tier2_ref': '5.6.1.2',
        'commentary_ref': 'A.4.1.2',
        'diaphragm_type': 'Flexible'
    },
    'STRAIGHT_SHEATHING': {
        'seismicity': ['High'],
        'statement': 'All straight-sheathed diaphragms have aspect ratios less than 2-to-1 in the direction being considered.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.1',
        'diaphragm_type': 'Flexible'
    },
    'SPANS': {
        'seismicity': ['High'],
        'statement': 'All wood diaphragms with spans greater than 24 ft (7.3 m) consist of wood structural panels or diagonal sheathing.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.2',
        'diaphragm_type': 'Flexible'
    },
    'DIAGONALLY_SHEATHED': {
        'seismicity': ['High'],
        'statement': 'All diagonally sheathed or unblocked wood structural panel diaphragms have horizontal spans less than 40 ft (12.2 m) and aspect ratios less than or equal to 4-to-1.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.3',
        'diaphragm_type': 'Flexible'
    },
    'OTHER_DIAPHRAGMS': {
        'seismicity': ['High'],
        'statement': 'Diaphragms do not consist of a system other than wood, metal deck, concrete, or horizontal bracing.',
        'tier2_ref': '5.6.5',
        'commentary_ref': 'A.4.7.1',
        'diaphragm_type': 'Flexible'
    },
    # Connections
    'UPLIFT_AT_PILE_CAPS': {
        'seismicity': ['High'],
        'statement': 'Pile caps have top reinforcement, and piles are anchored to the pile caps.',
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    }
}


# Table 17-25: Immediate Occupancy Structural Checklist for Building Types C2 and C2a
C2_IO_CHECKLIST = {
    # Very Low Seismicity (and also apply to Low/Moderate/High unless superseded)
    'COMPLETE_FRAMES': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Steel or concrete frames classified as secondary components form a complete vertical-load-carrying system.',
        'tier2_ref': '5.5.2.5.1',
        'commentary_ref': 'A.3.1.6.1'
    },
    'REDUNDANCY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The number of lines of shear walls in each principal direction is greater than or equal to 2.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.2.1.1'
    },
    'SHEAR_STRESS_CHECK': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': "The shear stress in the concrete shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than the greater of 100 lb/in² (0.69 MPa) or 2√f\'c.",
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.2.1'
    },
    'REINFORCING_STEEL': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The ratio of reinforcing steel area to gross concrete area is not less than 0.0012 in the vertical direction and 0.0020 in the horizontal direction. The spacing of reinforcing steel is equal to or less than 18 in. (457 mm).',
        'tier2_ref': '5.5.3.1.3',
        'commentary_ref': 'A.3.2.2.2'
    },
    'WALL_ANCHORAGE_FLEXIBLE': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Exterior concrete or masonry walls that are dependent on flexible diaphragms for lateral support are anchored for out-of-plane forces at each diaphragm level with steel anchors, reinforcing dowels, or straps that are developed into the diaphragm. Connections have strength to resist the connection force calculated in the Quick Check procedure of Section 4.4.3.7.',
        'tier2_ref': '5.7.1.1',
        'commentary_ref': 'A.5.1.1',
        'diaphragm_type': 'Flexible'
    },
    'TRANSFER_TO_SHEAR_WALLS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Diaphragms are connected for transfer of loads to the shear walls, and the connections are able to develop the lesser of the shear strength of the walls or diaphragms.',
        'tier2_ref': '5.7.2',
        'commentary_ref': 'A.5.2.1'
    },
    'FOUNDATION_DOWELS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Wall reinforcement is doweled into the foundation, and the dowels are able to develop the lesser of the strength of the walls or the uplift capacity of the foundation.',
        'tier2_ref': '5.7.3.4',
        'commentary_ref': 'A.5.3.5'
    },
    # Foundation system
    'DEEP_FOUNDATIONS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Piles and piers are capable of transferring the lateral forces between the structure and the soil.',
        'tier2_ref': '',
        'commentary_ref': 'A.6.2.3'
    },
    'SLOPING_SITES': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The difference in foundation embedment depth from one side of the building to another does not exceed one story.',
        'tier2_ref': '',
        'commentary_ref': 'A.6.2.4'
    },
    
    # Low, Moderate, and High Seismicity (additional items)
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components and are compliant with the following items: COLUMN-BAR SPLICES, BEAM-BAR SPLICES, COLUMN-TIE SPACING, STIRRUP SPACING, and STIRRUP AND TIE HOOKS.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Flat slabs or plates not part of seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'COUPLING_BEAMS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The ends of both walls to which the coupling beam is attached are supported at each end to resist vertical loads caused by overturning. Coupling beams have the capacity in shear to develop the uplift capacity of the adjacent wall.',
        'tier2_ref': '5.5.3.2.1',
        'commentary_ref': 'A.3.2.2.3'
    },
    'OVERTURNING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All shear walls have aspect ratios less than 4-to-1. Wall piers need not be considered.',
        'tier2_ref': '5.5.3.1.4',
        'commentary_ref': 'A.3.2.2.4'
    },
    'CONFINEMENT_REINFORCING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'For shear walls with aspect ratios greater than 2-to-1, the boundary elements are confined with spirals or ties with spacing less than 8db.',
        'tier2_ref': '5.5.3.2.2',
        'commentary_ref': 'A.3.2.2.5'
    },
    'WALL_REINFORCING_AT_OPENINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is added trim reinforcement around all wall openings with a dimension greater than three times the thickness of the wall.',
        'tier2_ref': '5.5.3.1.5',
        'commentary_ref': 'A.3.2.2.6'
    },
    'WALL_THICKNESS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Thicknesses of bearing walls are not less than 1/25 the unsupported height or length, whichever is shorter, nor less than 4 in. (101 mm).',
        'tier2_ref': '5.5.3.1.2',
        'commentary_ref': 'A.3.2.2.7'
    },
    # Diaphragms (Stiff or Flexible)
    'DIAPHRAGM_CONTINUITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The diaphragms are not composed of split-level floors and do not have expansion joints.',
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'OPENINGS_AT_SHEAR_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragm openings immediately adjacent to the shear walls are less than 15% of the wall length.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.4'
    },
    'PLAN_IRREGULARITIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is tensile capacity to develop the strength of the diaphragm at reentrant corners or other locations of plan irregularities.',
        'tier2_ref': '5.6.1.4',
        'commentary_ref': 'A.4.1.7'
    },
    'DIAPHRAGM_REINFORCEMENT_AT_OPENINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is reinforcing around all diaphragm openings larger than 50% of the building width in either major plan dimension.',
        'tier2_ref': '5.6.1.5',
        'commentary_ref': 'A.4.1.8'
    },
    # Flexible Diaphragms
    'CROSS_TIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There are continuous cross ties between diaphragm chords.',
        'tier2_ref': '5.6.1.2',
        'commentary_ref': 'A.4.1.2',
        'diaphragm_type': 'Flexible'
    },
    'STRAIGHT_SHEATHING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All straight-sheathed diaphragms have aspect ratios less than 1-to-1 in the direction being considered.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.1',
        'diaphragm_type': 'Flexible'
    },
    'SPANS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All wood diaphragms with spans greater than 12 ft (3.6 m) consist of wood structural panels or diagonal sheathing.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.2',
        'diaphragm_type': 'Flexible'
    },
    'DIAGONALLY_SHEATHED': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All diagonally sheathed or unblocked wood structural panel diaphragms have horizontal spans less than 30 ft (9.2 m) and aspect ratios less than or equal to 3-to-1.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.3',
        'diaphragm_type': 'Flexible'
    },
    'NONCONCRETE_FILLED': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Untopped metal deck diaphragms or metal deck diaphragms with fill other than concrete consist of horizontal spans of less than 40 ft (12.2 m) and have aspect ratios less than 4-to-1.',
        'tier2_ref': '5.6.3',
        'commentary_ref': 'A.4.3.1',
        'diaphragm_type': 'Flexible'
    },
    'OTHER_DIAPHRAGMS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragms do not consist of a system other than wood, metal deck, concrete, or horizontal bracing.',
        'tier2_ref': '5.6.5',
        'commentary_ref': 'A.4.7.1',
        'diaphragm_type': 'Flexible'
    },
    # Connections
    'UPLIFT_AT_PILE_CAPS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Pile caps have top reinforcement, and piles are anchored to the pile caps; the pile cap reinforcement and pile anchorage are able to develop the tensile capacity of the piles.',
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    }
}


# ============================================================================
# C3: CONCRETE FRAMES WITH INFILL MASONRY SHEAR WALLS
# C3A: CONCRETE FRAMES WITH INFILL MASONRY SHEAR WALLS AND FLEXIBLE DIAPHRAGMS
# ============================================================================

# Table 17-26: Collapse Prevention Structural Checklist for Building Types C3 and C3a
C3_CP_CHECKLIST = {
	# Low & Moderate Seismicity
    'REDUNDANCY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The number of lines of shear walls in each principal direction is greater than or equal to 2. ',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.2.1.1'
    },
    'SHEAR_STRESS_CHECK_1': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The shear stress in the reinforced masonry shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than the greater of 70 lb/in² (0.48 MPa)",
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.4.1'
    },
    'SHEAR_STRESS_CHECK_2': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': "The shear stress in the unreinforced masonry shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than 30 lb/in.2 (0.21 MPa) for clay units and 70 lb/in.2 (0.48 MPa) for concrete units. Bays with openings greater than 25% of the wall area shall not be included in Aw of Eq. (4-8).",
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.5.1'
    },
    'INFILL_WALL_CONNECTIONS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Masonry is in full contact with frame.',
        'tier2_ref': '5.5.3.5.1, 5.5.3.5.3',
        'commentary_ref': 'A.3.2.6.1'
    },
	# Connections
    'TRANSFER_TO_SHEAR_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragms are connected for transfer of loads to the shear walls, and the connections are able to develop the lesser of the shear strength of the walls or diaphragms.',
        'tier2_ref': '5.7.2',
        'commentary_ref': 'A.5.2.1'
    },
    'CONCRETE_COLUMNS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All concrete columns are doweled into the foundation with a minimum of four bars, and the dowels are able to develop the tensile capacity of reinforcement in columns of the seismic-force-resisting system.',
        'tier2_ref': '5.7.3.1',
        'commentary_ref': 'A.5.3.2'
    },
    # High Seismicity
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components and are Compliant with the following items: COLUMN-BAR SPLICES, BEAM-BAR SPLICES, COLUMN-TIE SPACING, STIRRUP SPACING, and STIRRUP AND TIE HOOKS.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['High'],
        'statement': 'Flat slabs or plates not part of the seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'PROPORTIONS': {
        'seismicity': ['High'],
        'statement': 'The height-to-thickness ratio of the unreinforced infill walls at each story is less than 13.',
        'tier2_ref': '5.5.3.1.2',
        'commentary_ref': 'A.3.2.6.2'
    },
    'CAVITY_WALLS': {
        'seismicity': ['High'],
        'statement': 'The infill walls are not of cavity construction.',
        'tier2_ref': '5.5.3.5.2',
        'commentary_ref': 'A.3.2.6.3'
    },
    'INFILL_WALLS': {
        'seismicity': ['High'],
        'statement': 'The infill walls are continuous to the soffits of the frame beams and to the columns to either side.',
        'tier2_ref': '5.5.3.5.3',
        'commentary_ref': 'A.3.2.6.4'
    },
	# Diaphragms (Stiff or Flexible)
    'DIAPHRAGM_CONTINUITY': {
        'seismicity': ['High'],
        'statement': 'The diaphragms are not composed of split-level floors and do not have expansion joints.',
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'OPENINGS_AT_SHEAR_WALLS': {
        'seismicity': ['High'],
        'statement': 'Diaphragm openings immediately adjacent to the shear walls are less than 15% of the wall length.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.4'
    },
    'OPENINGS_AT_EXTERIOR_WALLS': {
        'seismicity': ['High'],
        'statement': 'Diaphragm openings immediately adjacent to exterior masonry shear walls are not greater than 4 ft long.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.6'
    },
    # Flexible Diaphragms	
    'CROSS_TIES': {
        'seismicity': ['High'],
        'statement': 'There are continuous cross ties between diaphragm chords.',
        'tier2_ref': '5.6.1.2',
        'commentary_ref': 'A.4.1.2',
        'diaphragm_type': 'Flexible'
    },
    'STRAIGHT_SHEATHING': {
        'seismicity': ['High'],
        'statement': 'All straight-sheathed diaphragms have aspect ratios less than 1-to-1 in the direction being considered.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.1',
        'diaphragm_type': 'Flexible'
    },
    'SPANS': {
        'seismicity': ['High'],
        'statement': 'All wood diaphragms with spans greater than 12 ft (3.6 m) consist of wood structural panels or diagonal sheathing.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.2',
        'diaphragm_type': 'Flexible'
    },
    'DIAGONALLY_SHEATHED': {
        'seismicity': ['High'],
        'statement': 'All diagonally sheathed or unblocked wood structural panel diaphragms have horizontal spans less than 30 ft (9.2 m) and aspect ratios less than or equal to 3-to-1.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.3',
        'diaphragm_type': 'Flexible'
    },
    'OTHER_DIAPHRAGMS': {
        'seismicity': ['High'],
        'statement': 'Diaphragms do not consist of a system other than wood, metal deck, concrete, or horizontal bracing.',
        'tier2_ref': '5.6.5',
        'commentary_ref': 'A.4.7.1',
        'diaphragm_type': 'Flexible'
    },
    # Connections	
    'UPLIFT_AT_PILE_CAPS': {
        'seismicity': ['High'],
        'statement': 'Pile caps have top reinforcement, and piles are anchored to the pile caps; the pile cap reinforcement and pile anchorage are able to develop the tensile capacity of the piles.',
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    },
    'STIFFNESS_OF_WALL_ANCHORS': {
        'seismicity': ['High'],
        'statement': 'Anchors of concrete or masonry walls to wood structural elements are installed taut and are stiff enough to limit the relative movement between the wall and the diaphragm to no greater than 1/8 in. before engagement of the anchors.',
        'tier2_ref': '5.7.1.2',
        'commentary_ref': 'A.5.1.4'
    },
}


# Table 17-27: Immediate Occupancy Structural Checklist for Building Types C3 and C3a
C3_IO_CHECKLIST = {
    # Very Low Seismicity (and also apply to Low/Moderate/High unless superseded)
    'REDUNDANCY': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The number of lines of shear walls in each principal direction is greater than or equal to 2.',
        'tier2_ref': '5.5.1.1',
        'commentary_ref': 'A.3.2.1.1'
    },
    'SHEAR_STRESS_CHECK_REINFORCED': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The shear stress in the reinforced masonry shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than 70 lb/in² (0.48 MPa).',
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.4.1'
    },
    'SHEAR_STRESS_CHECK_UNREINFORCED': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'The shear stress in the unreinforced masonry shear walls, calculated using the Quick Check procedure of Section 4.4.3.3, is less than 30 lb/in² (0.21 MPa) for clay units and 70 lb/in² (0.48 MPa) for concrete units. Bays with openings greater than 25% of the wall area shall not be included in Aw of Eq. (4-8).',
        'tier2_ref': '5.5.3.1.1',
        'commentary_ref': 'A.3.2.5.1'
    },
    'INFILL_WALL_CONNECTIONS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Masonry is in full contact with frame.',
        'tier2_ref': '5.5.3.5.1, 5.5.3.5.3',
        'commentary_ref': 'A.3.2.6.1'
    },
    # Connections 
    'TRANSFER_TO_SHEAR_WALLS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'Diaphragms are connected for transfer of loads to the shear walls, and the connections are able to develop the lesser of the shear strength of the walls or diaphragms.',
        'tier2_ref': '5.7.2',
        'commentary_ref': 'A.5.2.1'
    },
    'CONCRETE_COLUMNS': {
        'seismicity': ['Very Low', 'Low', 'Moderate', 'High'],
        'statement': 'All concrete columns are doweled into the foundation with a minimum of four bars, and the dowels are able to develop the tensile capacity of reinforcement in columns of the seismic-force-resisting system.',
        'tier2_ref': '5.7.3.1',
        'commentary_ref': 'A.5.3.2'
    },
    
    # Low and Moderate Seismicity (additional items)
    # Seismic Force-resisting system
    'DEFLECTION_COMPATIBILITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Secondary components have the shear capacity to develop the flexural strength of the components and are Compliant with the following items: COLUMN-BAR SPLICES, BEAM-BAR SPLICES, COLUMN-TIE SPACING, STIRRUP SPACING, and STIRRUP AND TIE HOOKS.',
        'tier2_ref': '5.5.2.5.2',
        'commentary_ref': 'A.3.1.6.2'
    },
    'FLAT_SLABS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Flat slabs or plates not part of the seismic-force-resisting system have continuous bottom steel through the column joints.',
        'tier2_ref': '5.5.2.5.3',
        'commentary_ref': 'A.3.1.6.3'
    },
    'REINFORCING_AT_WALL_OPENINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All wall openings that interrupt rebar have trim reinforcing on all sides.',
        'tier2_ref': '5.5.3.1.5',
        'commentary_ref': 'A.3.2.4.3'
    },
    'PROPORTIONS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The height-to-thickness ratio of the unreinforced infill walls at each story is less than 13.',
        'tier2_ref': '5.5.3.1.2',
        'commentary_ref': 'A.3.2.6.2'
    },
    'CAVITY_WALLS': {
        'seismicity': ['Low', 'Moderate'],
        'statement': 'The infill walls are not of cavity construction.',
        'tier2_ref': '5.5.3.5.2',
        'commentary_ref': 'A.3.2.6.3'
    },
    'INFILL_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The infill walls are continuous to the soffits of the frame beams and to the columns to either side.',
        'tier2_ref': '5.5.3.5.3',
        'commentary_ref': 'A.3.2.6.4'
    },
    'DIAPHRAGM_CONTINUITY': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'The diaphragms are not composed of split-level floors and do not have expansion joints.',
        'tier2_ref': '5.6.1.1',
        'commentary_ref': 'A.4.1.1'
    },
    'OPENINGS_AT_SHEAR_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragm openings immediately adjacent to the shear walls are less than 15% of the wall length.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.4'
    },
    'OPENINGS_AT_EXTERIOR_WALLS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragm openings immediately adjacent to exterior masonry shear walls are not greater than 4 ft long.',
        'tier2_ref': '5.6.1.3',
        'commentary_ref': 'A.4.1.6'
    },
    'PLAN_IRREGULARITIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is tensile capacity to develop the strength of the diaphragm at reentrant corners or other locations of plan irregularities.',
        'tier2_ref': '5.6.1.4',
        'commentary_ref': 'A.4.1.7'
    },
    'DIAPHRAGM_REINFORCEMENT_AT_OPENINGS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There is reinforcing around all diaphragm openings larger than 50% of the building width in either major plan dimension.',
        'tier2_ref': '5.6.1.5',
        'commentary_ref': 'A.4.1.8'
    },
    # Flexible Diaphragms
    'CROSS_TIES': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'There are continuous cross ties between diaphragm chords.',
        'tier2_ref': '5.6.1.2',
        'commentary_ref': 'A.4.1.2',
        'diaphragm_type': 'Flexible'
    },
    'STRAIGHT_SHEATHING': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All straight-sheathed diaphragms have aspect ratios less than 1-to-1 in the direction being considered.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.1',
        'diaphragm_type': 'Flexible'
    },
    'SPANS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All wood diaphragms with spans greater than 12 ft (3.6 m) consist of wood structural panels or diagonal sheathing.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.2',
        'diaphragm_type': 'Flexible'
    },
    'DIAGONALLY_SHEATHED': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'All diagonally sheathed or unblocked wood structural panel diaphragms have horizontal spans less than 30 ft (9.2 m) and aspect ratios less than or equal to 3-to-1.',
        'tier2_ref': '5.6.2',
        'commentary_ref': 'A.4.2.3',
        'diaphragm_type': 'Flexible'
    },
    'NONCONCRETE_FILLED': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Untopped metal deck diaphragms or metal deck diaphragms with fill other than concrete consist of horizontal spans of less than 40 ft (12.2 m) and have aspect ratios less than 4-to-1.',
        'tier2_ref': '5.6.3',
        'commentary_ref': 'A.4.3.1',
        'diaphragm_type': 'Flexible'
    },
    'OTHER_DIAPHRAGMS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Diaphragms do not consist of a system other than wood, metal deck, concrete, or horizontal bracing.',
        'tier2_ref': '5.6.5',
        'commentary_ref': 'A.4.7.1',
        'diaphragm_type': 'Flexible'
    },
    # Connections
    'UPLIFT_AT_PILE_CAPS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Pile caps have top reinforcement, and piles are anchored to the pile caps; the pile cap reinforcement and pile anchorage are able to develop the tensile capacity of the piles.',
        'tier2_ref': '5.7.3.5',
        'commentary_ref': 'A.5.3.8'
    },
    'STIFFNESS_OF_WALL_ANCHORS': {
        'seismicity': ['Low', 'Moderate', 'High'],
        'statement': 'Anchors of concrete or masonry walls to wood structural elements are installed taut and are stiff enough to limit the relative movement between the wall and the diaphragm to no greater than 1/8 in. before engagement of the anchors.',
        'tier2_ref': '5.7.1.2',
        'commentary_ref': 'A.5.1.4'
    },
    
    # High Seismicity (additional items)
    'PROPORTIONS_HIGH': {
        'seismicity': ['High'],
        'statement': 'The height-to-thickness ratio of the unreinforced infill walls at each story is less than 8.',
        'tier2_ref': '5.5.3.1.2',
        'commentary_ref': 'A.3.2.6.2'
    }
}


# ============================================================================
# Master dictionary to access checklists by building type and performance level
# ============================================================================

TIER1_CHECKLISTS = {
    'all': {
        'CP': CP_BASIC_CONFIGURATION_CHECKLIST,
        'IO': IO_BASIC_CONFIGURATION_CHECKLIST,
    },
    'C1': {
        'CP': C1_CP_CHECKLIST,
        'IO': C1_IO_CHECKLIST,
    },
    'C2': {
        'CP': C2_CP_CHECKLIST,
        'IO': C2_IO_CHECKLIST,
    },
    'C2a': {
        'CP': C2_CP_CHECKLIST,
        'IO': C2_IO_CHECKLIST,
    },
    'C3': {
        'CP': C3_CP_CHECKLIST,
        'IO': C3_IO_CHECKLIST,
    },
    'C3a': {
        'CP': C3_CP_CHECKLIST,
        'IO': C3_IO_CHECKLIST,
    }
}


import csv
from typing import Dict, List, Any, Optional
import os


def generate_checklist_spreadsheet(
    seismicity: str,
    performance_level: str,
    output_file: Optional[str] = None,
    building_types: Optional[List[str]] = None,
    diaphragm_type: Optional[str] = 'Stiff',
) -> List[Dict[str, Any]]:
    """
    Generate a spreadsheet of checklist items filtered by seismicity and performance level.
    
    Args:
        seismicity: Seismicity level ('Very Low', 'Low', 'Moderate', 'High')
        performance_level: Performance level ('CP' or 'IO')
        output_file: Optional path to output CSV file. If None, returns data only.
        building_types: Optional list of building types to include (e.g., ['C1', 'C2']).
                       If None, includes all building types in TIER1_CHECKLISTS.
        diaphragm_type: Optional diaphragm type to include (e.g., 'Flexible' or 'Stiff').
    
    Returns:
        List of dictionaries, each representing a row in the spreadsheet.
    
    The spreadsheet will have the following columns:
        'Checklist', 'key', 'Tier2_ref', 'commentary_ref', 'notes', 'statement'
    """
    
    rows = []
    
    # Determine which building types to process
    if building_types is None:
        # Process all building types except 'all' (basic config)
        types_to_process = [bt for bt in TIER1_CHECKLISTS.keys() if bt != 'all']
        # Always include basic config
        include_basic = True
    else:
        types_to_process = [bt for bt in building_types if bt in TIER1_CHECKLISTS and bt != 'all']
        include_basic = 'all' in building_types if building_types else False
    
    # Process basic configuration checklists (if requested)
    if include_basic and 'all' in TIER1_CHECKLISTS:
        basic_checklist = TIER1_CHECKLISTS['all'].get(performance_level, {})
        for key, item in basic_checklist.items():
            # Check if this item applies to the given seismicity
            item_seismicity = item.get('seismicity', [])
            # if building seismicity matches the item list or if there is no item list
            if (seismicity in item_seismicity) or (not item_seismicity):
                rows.append({
                    'Checklist': f'all_{performance_level}',
                    'Key': key,
                    'Tier2_ref': item.get('tier2_ref', ''),
                    'Commentary_ref': item.get('commentary_ref', ''),
                    'Notes': '',
                    'Statement': item.get('statement', '')
                })
    
    # Process building-type-specific checklists
    for building_type in types_to_process:
        checklist_dict = TIER1_CHECKLISTS.get(building_type, {})
        performance_checklist = checklist_dict.get(performance_level, {})
        
        for key, item in performance_checklist.items():
            # Check if this item applies to the given seismicity
            item_seismicity = item.get('seismicity', [])
            item_diaphragm_type = item.get('diaphragm_type', diaphragm_type)
            if ((seismicity in item_seismicity) or (not item_seismicity)) and \
                (diaphragm_type == item_diaphragm_type):
                rows.append({
                    'Checklist': f'{building_type}_{performance_level}',
                    'Key': key,
                    'Tier2_ref': item.get('tier2_ref', ''),
                    'Commentary_ref': item.get('commentary_ref', ''),
                    'Notes': '',
                    'Statement': item.get('statement', '')
                })
    
    # Sort rows by Checklist, then by key
    # rows.sort(key=lambda x: (x['Checklist'], x['Key']))
    
    # Write to CSV if output_file is provided
    if output_file:
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Checklist', 'Key', 'Tier2_ref', 'Commentary_ref', 'Notes', 'Statement']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for row in rows:
                # Clean up any None values
                clean_row = {k: (v if v is not None else '') for k, v in row.items()}
                writer.writerow(clean_row)
        
        print(f"Spreadsheet written to: {output_file}")
        print(f"Total rows: {len(rows)}")
    
    return rows


def print_checklist_summary(rows: List[Dict[str, Any]]) -> None:
    """
    Print a summary of the generated checklist rows.
    
    Args:
        rows: List of rows returned by generate_checklist_spreadsheet
    """
    from collections import Counter
    
    checklist_counts = Counter(row['Checklist'] for row in rows)
    
    print("\n" + "="*60)
    print("CHECKLIST SUMMARY")
    print("="*60)
    for checklist, count in sorted(checklist_counts.items()):
        print(f"{checklist}: {count} items")
    print("="*60)
    print(f"TOTAL: {len(rows)} items")
    print("="*60)


def generate_checklist_excel(
    seismicity: str,
    performance_level: str,
    output_file: str,
    building_types: Optional[list[str]] = None
) -> None:
    """
    Generate an Excel spreadsheet of checklist items.
    
    Requires openpyxl: pip install openpyxl
    """
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    except ImportError:
        print("openpyxl not installed. Please install it: pip install openpyxl")
        return
    
    # Get the data
    rows = generate_checklist_spreadsheet(seismicity, performance_level, None, building_types)
    
    # Create workbook
    wb = Workbook()
    ws = wb.active
    ws.title = f"Checklists_{seismicity}_{performance_level}"
    
    # Define headers
    headers = ['Checklist', 'Key', 'Tier2 Ref', 'Commentary Ref', 'Notes', 'Statement']
    ws.append(headers)
    
    # Style headers
    header_font = Font(bold=True, color="FFFFFF")
    header_fill = PatternFill(start_color="366092", end_color="366092", fill_type="solid")
    header_alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    
    for col in range(1, 7):
        cell = ws.cell(row=1, column=col)
        cell.font = header_font
        cell.fill = header_fill
        cell.alignment = header_alignment
    
    # Add data
    thin_border = Border(
        left=Side(style='thin'),
        right=Side(style='thin'),
        top=Side(style='thin'),
        bottom=Side(style='thin')
    )
    
    for row in rows:
        ws.append([
            row['Checklist'],
            row['key'],
            row['Tier2_ref'],
            row['commentary_ref'],
            row['notes'],
            row['statement']
        ])
        
        # Apply border to the new row
        for col in range(1, 7):
            ws.cell(row=ws.max_row, column=col).border = thin_border
    
    # Adjust column widths
    column_widths = [15, 20, 12, 15, 10, 80]
    for i, width in enumerate(column_widths, 1):
        ws.column_dimensions[chr(64 + i)].width = width
    
    # Freeze header row
    ws.freeze_panes = 'A2'
    
    # Save
    wb.save(output_file)
    print(f"Excel spreadsheet written to: {output_file}")
    print(f"Total rows: {len(rows)}")


# ------------------------------------------------------------------
# Sample building data (SI units: m, kN, MPa)
# ------------------------------------------------------------------

building_C1_data = {
    'Ss': 0.375,
    'S1': 0.15,
    'site_class': 'B',
    # 'seismic_design_category': 'D',
    # 'R': 8.0,             # Response modification coefficient (optional)
    # 'I_e': 1.0,           # Importance factor
    # 'use_code_period': True,  # Use approximate period
    'structure_type': 'C1',
    
    'performance_level': 'CP',
    # 'seismicity_level': 'Moderate',
    # 'building_height_m': 12.0,
    # 'base_shear_kN': 1500.0,          # pseudo‑lateral force (Section 4.4.2)
    'storey_elevations_m': [12.0, 8.0, 4.0],   # top to bottom
    'storey_weights_kN': [1000.0, 1200.0, 1300.0],
    # Direction‑specific data (both X and Y assumed identical for simplicity)
    'directions': {
        'X': {
            'diaphragm': 'rigid',
            # 'story_shear_kN': [500.0, 1000.0, 1500.0],   # top to bottom
            'total_column_area_m2': 0.8,          # total column cross‑section area
            'num_columns': 12,
            'num_frames': 3,
            'frame_length_m': 18.0,
            'end_column_area_m2': 0.08,           # area of a single end column
            'f_c_MPa': 20.0,
            'joints': {                            # optional, for joint check
                'joint_width_m': 0.4,
                'column_depth_m': 0.5,
                'beam_depth_m': 0.4
            }
        },
        'Y': {
            'diaphragm': 'rigid',
            # 'story_shear_kN': [500.0, 1000.0, 1500.0],
            'total_column_area_m2': 0.8,
            'num_columns': 12,
            'num_frames': 3,
            'frame_length_m': 18.0,
            'end_column_area_m2': 0.08,
            'f_c_MPa': 25.0,
            'joints': {
                'joint_width_m': 0.4,
                'column_depth_m': 0.5,
                'beam_depth_m': 0.4
            }
        }
    },
    # Configuration checks (true = compliant, false = non‑compliant)
    'config_checks': {
        'load_path': True,
        'adjacent_buildings': True,
        'meandering': True,
        'vertical_irregular': False,   # example: soft storey present
        'plan_irregular': True,
        'diaphragm_continuity': True
    }
}

building_C2_data = {
  "structure_type": "C2",
  "performance_level": "CP",
  "seismicity_level": "Moderate",
  "building_name": "Shear Wall Building",
#   "building_height_m": 20.0,
  "storey_elevations_m": [24.175, 20.075, 15.975, 11.550, 7.125],
  "storey_weights_kN": [800 + 7467.0, 10476.0, 16317.0, 17699.0, 0.0],
  "Ss": 0.375,
  "S1": 0.15,
  "site_class": "B",
  "risk_category": "IV",
  "seismic_design_level": "special",
  "I_e": 1.0,
  "use_code_period": True,
  
  "concrete_properties": {
    "f_c_MPa": 20.0,
    "E_c_MPa": 26000,
    "rebar_fy_MPa": 460,
    "rebar_fu_MPa": 620
  },
  
  "wall_system": {
    "typical_thickness_m": 0.15,
    "boundary_element_confined": False,
    "confining_tie_spacing_mm": 1000,
    "vertical_reinf_ratio": 0.0040,
    "horizontal_reinf_ratio": 0.0020,
    "vertical_bar_size_mm": 12,
    "horizontal_bar_size_mm": 12,
    "vertical_bar_spacing_mm": 300,
    "horizontal_bar_spacing_mm": 300,
    "trim_reinforcement_at_openings": True
  },
  
  "frame_system": {
    "secondary_frames_present": True,
    "beam_info": {
      "typical_depth_m": 0.75,
      "typical_width_m": 0.4,
      "longitudinal_bars_top": 5,
      "longitudinal_bars_bottom": 3,
      "bar_size_mm": 32,
      "stirrup_size_mm": 12,
      "stirrup_spacing_mm": 125
    },
    "column_info": {
      "typical_dimension_m": 0.4,
      "longitudinal_bars": 4,
      "bar_size_mm": 25,
      "tie_size_mm": 8,
      "tie_spacing_mm": 300,
      "tie_angle_degrees": 90
    }
  },
  
  "directions": {
    "X": {
      "num_frames": 4,
      "num_columns": 12,
      "frame_length_m": 6.0,
      "total_column_area_m2": 13*5*0.4*0.4,
      "end_column_area_m2": 0.16,
      "f_c_MPa": 20.0,
      "story_shear_kN": [400.0, 800.0, 1200.0, 1600.0, 2000.0],
      "lines_of_shear_walls": 2,
      "wall_info": {
        "wall_1": {
          "length_m": 6.0,
          "thickness_m": 0.15,
          "area_m2": 0.6,
          "aspect_ratio": 4.2,
          "openings": []
        },
        "wall_2": {
          "length_m": 6.0,
          "thickness_m": 0.15,
          "area_m2": 0.6,
          "aspect_ratio": 4.2,
          "openings": []
        },
        "wall_3": {
          "length_m": 6.0,
          "thickness_m": 0.15,
          "area_m2": 0.6,
          "aspect_ratio": 4.2,
          "openings": []
        },
        "wall_4": {
          "length_m": 4.0,
          "thickness_m": 0.15,
          "area_m2": 0.6,
          "aspect_ratio": 4.2,
          "openings": []
        },
        "wall_5": {
          "length_m": 10.0,
          "thickness_m": 0.15,
          "area_m2": 1.5,
          "aspect_ratio": 1.75,
          "openings": []
        },
        "wall_6": {
          "length_m": 14.0,
          "thickness_m": 0.15,
          "area_m2": 2.1,
          "aspect_ratio": 1.22,
          "openings": []
        }
      },
      "total_wall_area_m2": 5.7,
      "wall_vertical_reinf_ratio": 0.0025,
      "wall_horizontal_reinf_ratio": 0.0030,
      "wall_vertical_bar_size_mm": 16,
      "wall_horizontal_bar_size_mm": 13,
      "wall_vertical_bar_spacing_mm": 300,
      "wall_horizontal_bar_spacing_mm": 300,
      "boundary_elements_confined": False,
      "confining_tie_spacing_mm": 300,
      "coupling_beams_present": False,
    #   "coupling_beam_details": {
    #     "span_m": 2.0,
    #     "depth_m": 0.6,
    #     "reinforcement_ratio": 0.015
    #   }
    },
    "Y": {
      "num_frames": 3,
      "num_columns": 12,
      "frame_length_m": 40.0,
      "total_column_area_m2": 1.92,
      "end_column_area_m2": 0.16,
      "f_c_MPa": 28.0,
      "story_shear_kN": [400.0, 800.0, 1200.0, 1600.0, 2000.0],
      "lines_of_shear_walls": 4,
      "wall_info": {
        "wall_1": {
          "length_m": 4.0,
          "thickness_m": 0.3,
          "area_m2": 2.1,
          "aspect_ratio": 2.86,
          "openings": []
        },
        "wall_2": {
          "length_m": 4.0,
          "thickness_m": 0.3,
          "area_m2": 2.1,
          "aspect_ratio": 2.86,
          "openings": []
        },
        "wall_3": {
          "length_m": 4.0,
          "thickness_m": 0.3,
          "area_m2": 1.5,
          "aspect_ratio": 4.0,
          "openings": []
        },
        "wall_4": {
          "length_m": 4.0,
          "thickness_m": 0.3,
          "area_m2": 1.5,
          "aspect_ratio": 4.0,
          "openings": []
        },
        "wall_5": {
          "length_m": 10.0,
          "thickness_m": 0.15,
          "area_m2": 1.5,
          "aspect_ratio": 1.75,
          "openings": []
        },
        "wall_6": {
          "length_m": 14.0,
          "thickness_m": 0.15,
          "area_m2": 2.1,
          "aspect_ratio": 1.22,
          "openings": []
        }   
      },
      "total_wall_area_m2": 7.2,
      "wall_vertical_reinf_ratio": 0.0025,
      "wall_horizontal_reinf_ratio": 0.0030,
      "wall_vertical_bar_size_mm": 16,
      "wall_horizontal_bar_size_mm": 13,
      "wall_vertical_bar_spacing_mm": 300,
      "wall_horizontal_bar_spacing_mm": 300,
      "boundary_elements_confined": False,
      "confining_tie_spacing_mm": 300,
      "coupling_beams_present": False
    }
  },
  
  "wall_information": [
    {
      "storey": 0,
      "orientation": "X",
      "thickness_m": 0.3,
      "length_m": 8.0,
      "height_m": 4.0,
      "type": "concrete",
      "aspect_ratio": 2.5,
      "vertical_reinf_ratio": 0.0025,
      "horizontal_reinf_ratio": 0.0030,
      "vertical_bar_size_mm": 16,
      "horizontal_bar_size_mm": 13,
      "vertical_bar_spacing_mm": 300,
      "horizontal_bar_spacing_mm": 300,
      "boundary_element_confined": True,
      "openings": []
    },
    {
      "storey": 0,
      "orientation": "X",
      "thickness_m": 0.3,
      "length_m": 6.0,
      "height_m": 4.0,
      "type": "concrete",
      "aspect_ratio": 3.33,
      "vertical_reinf_ratio": 0.0025,
      "horizontal_reinf_ratio": 0.0030,
      "vertical_bar_size_mm": 16,
      "horizontal_bar_size_mm": 13,
      "vertical_bar_spacing_mm": 300,
      "horizontal_bar_spacing_mm": 300,
      "boundary_element_confined": True,
      "openings": []
    },
    {
      "storey": 0,
      "orientation": "X",
      "thickness_m": 0.3,
      "length_m": 5.0,
      "height_m": 4.0,
      "type": "concrete",
      "aspect_ratio": 4.0,
      "vertical_reinf_ratio": 0.0025,
      "horizontal_reinf_ratio": 0.0030,
      "vertical_bar_size_mm": 16,
      "horizontal_bar_size_mm": 13,
      "vertical_bar_spacing_mm": 300,
      "horizontal_bar_spacing_mm": 300,
      "boundary_element_confined": False,
      "openings": []
    },
    {
      "storey": 0,
      "orientation": "Y",
      "thickness_m": 0.3,
      "length_m": 7.0,
      "height_m": 4.0,
      "type": "concrete",
      "aspect_ratio": 2.86,
      "vertical_reinf_ratio": 0.0025,
      "horizontal_reinf_ratio": 0.0030,
      "vertical_bar_size_mm": 16,
      "horizontal_bar_size_mm": 13,
      "vertical_bar_spacing_mm": 300,
      "horizontal_bar_spacing_mm": 300,
      "boundary_element_confined": True,
      "openings": []
    },
    {
      "storey": 0,
      "orientation": "Y",
      "thickness_m": 0.3,
      "length_m": 7.0,
      "height_m": 4.0,
      "type": "concrete",
      "aspect_ratio": 2.86,
      "vertical_reinf_ratio": 0.0025,
      "horizontal_reinf_ratio": 0.0030,
      "vertical_bar_size_mm": 16,
      "horizontal_bar_size_mm": 13,
      "vertical_bar_spacing_mm": 300,
      "horizontal_bar_spacing_mm": 300,
      "boundary_element_confined": True,
      "openings": []
    }
  ],
  
  "diaphragm_info": {
    "type": "concrete",
    "thickness_mm": 200,
    "continuous": True,
    "split_level": False,
    "expansion_joints": False,
    "openings": [
      {
        "location": "level_2",
        "size_percent": 30,
        "adjacent_to_wall": False,
        "wall_length_affected_m": 0
      }
    ],
    "chord_continuous": True,
    "span_m": 15.0,
    "aspect_ratio": 1.2,
    "cross_ties_present": True,
    "straight_sheathing": False,
    "diagonal_sheathing": False,
    "untopped_metal_deck": False
  },
  
  "foundation_info": {
    "type": "mat_foundation",
    "ties_present": True,
    "soil_class": "D",
    "liquefaction_potential": False,
    "slope_failure_potential": False,
    "fault_rupture_potential": False,
    "deep_foundations": False,
    "sloping_site": False,
    "foundation_embedment_diff_m": 0.0,
    "overturning_ratio": 1.2,
    "pile_caps": {
      "top_reinforcement": True,
      "piles_anchored": True,
      "tensile_capacity_developed": True
    },
    "wall_dowels": {
      "size_equal_to_wall_reinf": True,
      "spacing_equal_to_wall_reinf": True,
      "able_to_develop_strength": True
    }
  },
  
  "config_checks": {
    "load_path": True,
    "adjacent_buildings": True,
    "mezzanines": True,
    "weak_story": True,
    "soft_story": True,
    "vertical_irregularities": True,
    "geometry": True,
    "mass": True,
    "torsion": True,
    "liquefaction": True,
    "slope_failure": True,
    "surface_fault_rupture": True,
    "overturning": True,
    "ties_between_foundation_elements": True
  }
}


def main():
    """
    Example: Run Tier 1 checks on a sample C1 building.
    In practice, replace the hardcoded dictionary with data loaded from a JSON file.
    """
    building_data = building_C2_data

    # ------------------------------------------------------------------
    # Generate checklists
    # ------------------------------------------------------------------
    # Option 1: Generate for High seismicity, IO performance level
    # Option 2: Generate for Moderate seismicity, CP performance level   
    # Option 3: Generate for Low seismicity, CP performance level
    # Option 4: Generate for Very Low seismicity, CP performance level
    option = 2

    # Example 1: Generate for High seismicity, IO performance level
    if option == 1:
        rows_high_io = generate_checklist_spreadsheet(
            seismicity='High',
            performance_level='IO',
            output_file='./checklists_high_io.csv'
        )
        print_checklist_summary(rows_high_io)
    
    if option == 2:
        # Example 2: Generate for Moderate seismicity, CP performance level
        # Only for specific building types
        rows_mod_cp = generate_checklist_spreadsheet(
            seismicity='Moderate',
            performance_level='CP',
            output_file='./checklists_C2_cp.csv',
            building_types=['all', 'C2']
        )
        print_checklist_summary(rows_mod_cp)
    
    if option == 3:
        # Example 3: Generate for Low seismicity, CP performance level
        # Return data only (no file)
        rows_low_cp = generate_checklist_spreadsheet(
            seismicity='Low',
            performance_level='CP'
        )
        print(f"\nLow seismicity CP checklists: {len(rows_low_cp)} items")
    
    if option == 4:
        # Example 4: Generate for Very Low seismicity, CP performance level
        rows_vlow_cp = generate_checklist_spreadsheet(
            seismicity='Very Low',
            performance_level='CP'
        )
        print(f"Very Low seismicity CP checklists: {len(rows_vlow_cp)} items")

    # ------------------------------------------------------------------
    # Run the checks
    # ------------------------------------------------------------------
    checker = Tier1Checker(building_data, calculate_forces=True)
    
    # Print period and base shear results
    print("\n" + "="*70)
    print("SEISMIC LOAD CALCULATION RESULTS")
    print("="*70)
    period_results = checker.data['period_analysis']
    # print(f"Approximate period (Ta): {period_results['Ta']:.3f} s")
    # print(f"Period used: {period_results['T_used']:.3f} s ({period_results['period_source']})")
    print(f"Code Period: {period_results['T_code']:.3f} s")
    print(f"Base shear coefficient (Cs): {period_results['Cs']:.4f}")
    # print(f"Governing equation: {period_results['governing_eq']}")
    print(f"Total building weight: {checker.data['total_weight_kN']:,.0f} kN")
    print(f"Base shear: {checker.data['base_shear_kN']:,.0f} kN")
    
    # Print storey forces
    print("\n" + "-"*70)
    print("Storey Forces (top to bottom)")
    print("-"*70)
    print(f"{'Storey':<10} {'Elev (m)':<10} {'Weight (kN)':<12} {'Force (kN)':<12} {'Shear (kN)':<12} {'Moment (kNm)':<12}")
    print("-"*70)
    for i, (elev, wt, force, shear, moment) in enumerate(zip(
        checker.data['storey_elevations_m'],
        checker.data['storey_weights_kN'],
        checker.data['storey_forces_kN'],
        checker.data['storey_shears_kN'],
        checker.data['storey_moments_kNm']
    )):
        print(f"{i+1:<10} {elev:<10.2f} {wt:<12.0f} {force:<12.0f} {shear:<12.0f} {moment:<12.0f}")
    
    # Run the rest of the checks
    results = checker.run_checks()

    # ------------------------------------------------------------------
    # Print results in a simple table
    # ------------------------------------------------------------------
    print("\n" + "="*70)
    print("TIER 1 SCREENING RESULTS")
    print("="*70)
    print(f"Site class: {building_data['site_class']}")
    print(f"Structure type: {checker.struct_type}")
    print(f"Performance level: {checker.perf}")
    print(f"Seismicity level: {checker.seis}")   # or checker.data['seismicity_level']
    print(f"S_s = {checker.data['Ss']:.3f}g, S_1 = {checker.data['S1']:.3f}g")
    if 'S_DS' in checker.data:
        print(f"S_DS = {checker.data['S_DS']:.3f}g, S_D1 = {checker.data['S_D1']:.3f}g")
    # print(f"Structure type: {building_data['structure_type']}")
    # print(f"Performance level: {building_data['performance_level']}")
    # print(f"Seismicity: {building_data['seismicity_level']}")
    print("-"*70)

    for res in results:
        # Basic info
        check_name = res.get('check', 'Unknown check')
        direction = res.get('direction', '')
        if direction:
            check_name = f"{check_name} ({direction})"

        # Pass/fail indicator
        passed = res.get('pass')
        if passed is True:
            status = "✅ PASS"
        elif passed is False:
            status = "❌ FAIL"
        else:
            status = "⚠️  NOT CHECKED"

        # Demand/capacity if available
        dc = res.get('dc_ratio')
        if dc is not None:
            status += f"  D/C = {dc:.3f}"

        # Additional notes
        note = res.get('note', '')
        if note:
            status += f"  [{note}]"

        print(f"{check_name:<50} {status}")

    print("="*70)


if __name__ == "__main__":
    main()
