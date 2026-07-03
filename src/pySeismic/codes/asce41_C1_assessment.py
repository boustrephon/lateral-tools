"""
asce41_C1_assessment.py - ASCE41-17 assessment functions for concrete frames
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Union, Optional
from pathlib import Path

# Import your existing modules
from pySeismic.codes.ASCE41_17 import *
from pySeismic.utilities import *
from pySeismic.codes.ASCE7_16 import site_coeffs
from pySeismic.plotting import plot_storey_VM, tripartite_plot


class ASCE41Assessment:
    """Main class for ASCE41-17 seismic assessment"""
    
    def __init__(self, building_info: Dict):
        self.building_info = building_info
        self.hazard_params = None
        self.base_shear_results = None
        self.storey_results = None
        self.checks = []
        
    def summarize_building_info(self) -> Dict:
        """Create summary tables of building information"""
        info = self.building_info
        
        print("="*80)
        print("BUILDING INFORMATION SUMMARY")
        print("="*80)
        
        # Seismic Hazard Summary
        print("\n1. SEISMIC HAZARD")
        print("-"*40)
        print(f"Short Period Spectral Acceleration, S_s: {info.get('S_s', 'N/A')} g")
        print(f"1-second Spectral Acceleration, S_1: {info.get('S_1', 'N/A')} g")
        print(f"Site Class: {info.get('site_class', 'N/A')}")
        print(f"Soil Class: {info.get('soil_class', 'N/A')}")
        
        # Calculate site coefficients
        F_a, F_v = site_coeffs(info.get('site_class', 'D'), 
                              info.get('S_s', 1.0), 
                              info.get('S_1', 0.4))
        S_XS = F_a * info.get('S_s', 1.0)
        S_X1 = F_v * info.get('S_1', 0.4)
        
        print(f"Site Coefficient F_a: {F_a:.3f}")
        print(f"Site Coefficient F_v: {F_v:.3f}")
        print(f"MCE Spectral Acceleration, S_XS: {S_XS:.3f} g")
        print(f"MCE Spectral Acceleration, S_X1: {S_X1:.3f} g")
        
        # Seismicity Level
        seismicity = seismicity_level_get(S_XS, S_X1)
        print(f"Level of Seismicity: {seismicity}")
        
        self.hazard_params = {
            'S_XS': S_XS,
            'S_X1': S_X1,
            'seismicity_level': seismicity,
            'F_a': F_a,
            'F_v': F_v
        }
        
        return self.hazard_params
    
    def calculate_base_shear(self) -> Dict:
        """Calculate base shear according to ASCE41-17"""
        info = self.building_info
        
        # Calculate building period
        T_code = T_code_get(info.get('building_type', 'C1'), 
                           info.get('total_height', 30.0))
        
        S_s = info.get('S_s', 1.0)
        S_1 = info.get('S_1', 0.4)
        site_class = info.get('site_class', 'D')
        building_type = info.get('building_type', 'C1')
        
        # Calculate spectral acceleration at building period
        S_a = S_a_func(T_code, S_s, S_1, site_class)
        
        # Calculate total seismic weight
        total_weight = sum(info.get('storey_weights', [0]))
        
        # Calculate base shear (Tier 1 assessment)
        S_a_T1 = S_a_Tier_1_func(T_code, S_s, S_1, site_class)
        
        # Get modification factor C from Table 4-7
        n_storeys = info.get('num_storeys', 1)
        C_factor = table_4_7_dict_17_get(building_type, n_storeys)
        
        # Calculate base shear
        base_shear = S_a_T1 * total_weight * (C_factor if C_factor else 1.0)
        
        self.base_shear_results = {
            'T': T_code,
            'S_a': S_a,
            'S_a_T1': S_a_T1,
            'total_weight': total_weight,
            'C_factor': C_factor,
            'base_shear': base_shear
        }
        
        return self.base_shear_results
    
    def calculate_storey_forces(self, k: float = 1.0) -> Tuple:
        """Calculate storey force distribution"""
        info = self.building_info
        
        storey_elevs = info.get('storey_elevations', [])
        storey_weights = info.get('storey_weights', [])
        
        if not storey_elevs or not storey_weights:
            return None, None, None, None
        
        # Ensure data is from top to bottom
        if storey_elevs[0] > storey_elevs[-1]:  # Top to bottom
            pass
        else:  # Bottom to top, reverse
            storey_elevs = storey_elevs[::-1]
            storey_weights = storey_weights[::-1]
        
        V_base = self.base_shear_results['base_shear']
        
        # Calculate storey force coefficients
        coeffs = storey_force_coeffs(storey_elevs, storey_weights, k)
        
        # Calculate storey forces
        forces = storey_force_distribution(V_base, storey_elevs, storey_weights, k)
        
        # Calculate storey shears and moments
        storey_heights = [storey_elevs[i] - storey_elevs[i+1] for i in range(len(storey_elevs)-1)] + [0]
        shears, moments = storey_shear_moment(forces, storey_heights, ascending=False)
        
        # Create results DataFrame
        results = pd.DataFrame({
            'Storey': info.get('storey_names', [f'Level {i+1}' for i in range(len(storey_elevs))]),
            'Elevation (m)': storey_elevs,
            'Weight (kN)': storey_weights,
            'Force Coeff': [f'{c:.4f}' for c in coeffs],
            'Lateral Force (kN)': [f'{f:,.0f}' for f in forces],
            'Storey Shear (kN)': [f'{s:,.0f}' for s in shears],
            'Storey Moment (kNm)': [f'{m:,.0f}' for m in moments]
        })
        
        self.storey_results = results
        return results, forces, shears, moments
    
    def check_column_shear_stress(self, storey_shear: float, direction: str = 'X') -> Dict:
        """Check concrete frame column shear stress per ASCE41-17 Section 4.4.3.2"""
        info = self.building_info
        
        # Get required parameters
        V_j = storey_shear  # Storey shear (kN)
        A_c = info.get('total_column_area', 0)  # Total column area (m²)
        n_c = info.get(f'frame_column_number_{direction}', 0)
        n_f = info.get(f'frame_number_{direction}', 0)
        perf = info.get('performance_level', 'IO')
        f_c = info.get('f_c', 35)  # MPa
        
        # Convert to US customary units for ASCE41 functions
        V_j_kip = kN2kip(abs(V_j))
        A_c_sqin = sqm2sqin(A_c)
        
        # Calculate demand
        demand_psi = v_j_avg_rcmf(V_j_kip, A_c_sqin, n_c, n_f, perf)
        demand_MPa = ksi2MPa(demand_psi)
        
        # Calculate capacity (Table 17-23 for IO)
        if perf == 'IO':
            capacity_psi = max(2 * MPa2psi(f_c)**0.5, 100.0)
        elif perf == 'LS':
            capacity_psi = 4 * MPa2psi(f_c)**0.5
        elif perf == 'CP':
            capacity_psi = 6 * MPa2psi(f_c)**0.5
        else:
            capacity_psi = 2 * MPa2psi(f_c)**0.5
        
        capacity_MPa = psi2MPa(capacity_psi)
        
        # Calculate D/C ratio
        dc_ratio = demand_MPa / capacity_MPa if capacity_MPa > 0 else 999
        
        result = {
            'check': f'Column Shear Stress ({direction})',
            'section': '4.4.3.2',
            'equation': '4-7',
            'demand_MPa': demand_MPa,
            'capacity_MPa': capacity_MPa,
            'dc_ratio': dc_ratio,
            'status': 'PASS' if dc_ratio <= 1.0 else 'FAIL'
        }
        
        self.checks.append(result)
        return result
    
    def check_column_axial_stress(self, direction: str = 'X') -> Dict:
        """Check column axial stress due to overturning per ASCE41-17 Section 4.4.3.6"""
        info = self.building_info
        
        # Get required parameters
        V = self.base_shear_results['base_shear']  # Base shear (kN)
        h_n = info.get('total_height', 0)  # Building height (m)
        L = info.get(f'frame_length_{direction}', 0)  # Frame length (m)
        n_f = info.get(f'frame_number_{direction}', 0)
        A_col = info.get(f'end_column_area_{direction}', 0)  # End column area (m²)
        perf = info.get('performance_level', 'IO')
        f_c = info.get('f_c', 35)  # MPa
        
        # Convert to US customary units
        V_kip = kN2kip(abs(V))
        h_n_ft = m2ft(h_n)
        L_ft = m2ft(L)
        A_col_sqin = sqm2sqin(A_col)
        
        # Calculate demand
        demand_psi = p_ot(V_kip, h_n_ft, L_ft, n_f, A_col_sqin, perf)
        demand_MPa = ksi2MPa(demand_psi)
        
        # Calculate capacity (Table 17-23 for IO)
        if perf == 'IO':
            capacity_psi = 0.3 * MPa2psi(f_c)  # 30% of f'c
        elif perf == 'LS':
            capacity_psi = 0.5 * MPa2psi(f_c)  # 50% of f'c
        elif perf == 'CP':
            capacity_psi = 0.7 * MPa2psi(f_c)  # 70% of f'c
        else:
            capacity_psi = 0.3 * MPa2psi(f_c)
        
        capacity_MPa = psi2MPa(capacity_psi)
        
        # Calculate D/C ratio
        dc_ratio = demand_MPa / capacity_MPa if capacity_MPa > 0 else 999
        
        result = {
            'check': f'Column Axial Stress ({direction})',
            'section': '4.4.3.6',
            'equation': '4-11',
            'demand_MPa': demand_MPa,
            'capacity_MPa': capacity_MPa,
            'dc_ratio': dc_ratio,
            'status': 'PASS' if dc_ratio <= 1.0 else 'FAIL'
        }
        
        self.checks.append(result)
        return result
    
    def check_beam_column_joint_shear(self, direction: str = 'X') -> Dict:
        """Check beam-column joint shear per ASCE41-17 Section 10.4.2"""
        info = self.building_info
        
        # This is a simplified check - full implementation would require
        # detailed joint geometry and reinforcement details
        perf = info.get('performance_level', 'IO')
        f_c = info.get('f_c', 35)  # MPa
        
        # Simplified joint shear check
        # In practice, you'd need to calculate joint shear forces from frame analysis
        demand_MPa = 0.0  # Placeholder - would come from analysis
        capacity_MPa = 0.0  # Placeholder - would be calculated per ASCE41-17
        
        # Placeholder values for demonstration
        if perf == 'IO':
            capacity_MPa = 0.25 * f_c
        elif perf == 'LS':
            capacity_MPa = 0.5 * f_c
        elif perf == 'CP':
            capacity_MPa = 0.75 * f_c
        
        dc_ratio = demand_MPa / capacity_MPa if capacity_MPa > 0 else 0
        
        result = {
            'check': f'Beam-Column Joint Shear ({direction})',
            'section': '10.4.2',
            'equation': '10-7',
            'demand_MPa': demand_MPa,
            'capacity_MPa': capacity_MPa,
            'dc_ratio': dc_ratio,
            'status': 'PASS' if dc_ratio <= 1.0 else 'FAIL',
            'note': 'Simplified check - requires detailed joint analysis'
        }
        
        self.checks.append(result)
        return result
    
    def check_drift_limits(self) -> Dict:
        """Check story drift limits per ASCE41-17 Section 7.5"""
        info = self.building_info
        perf = info.get('performance_level', 'IO')
        
        # Drift limits from Table 7-10
        drift_limits = {
            'IO': 0.01,  # 1% of story height
            'LS': 0.02,  # 2% of story height
            'CP': 0.025  # 2.5% of story height
        }
        
        limit = drift_limits.get(perf, 0.02)
        
        # Placeholder - actual drift would come from analysis
        max_drift_ratio = 0.005  # Example value
        
        dc_ratio = max_drift_ratio / limit if limit > 0 else 0
        
        result = {
            'check': 'Story Drift Limits',
            'section': '7.5',
            'table': '7-10',
            'demand_ratio': max_drift_ratio,
            'limit_ratio': limit,
            'dc_ratio': dc_ratio,
            'status': 'PASS' if dc_ratio <= 1.0 else 'FAIL',
            'note': 'Placeholder - requires structural analysis for actual drifts'
        }
        
        self.checks.append(result)
        return result
    
    def check_tier1_limitations(self) -> bool:
        """Check Tier 1 procedure limitations per ASCE41-17 Table 3-4"""
        info = self.building_info
        
        building_type = info.get('building_type', 'C1')
        seismicity_level = self.hazard_params['seismicity_level']
        performance_level = info.get('performance_level', 'IO')
        num_storeys = info.get('num_storeys', 1)
        
        # Convert performance level to Table 3-4 format
        perf_table = {
            'IO': 'S-1',
            'LS': 'S-1',
            'CP': 'S-5'
        }.get(performance_level, 'S-1')
        
        # Get limitation from Table 3-4
        limit = table_3_4_get(building_type, seismicity_level, perf_table)
        
        if limit == 'NL':
            return True
        elif limit == 'NP':
            return False
        else:
            try:
                limit_int = int(limit)
                return num_storeys <= limit_int
            except:
                return None
    
    def generate_summary_report(self) -> str:
        """Generate comprehensive summary report"""
        report = []
        report.append("="*80)
        report.append("ASCE41-17 SEISMIC ASSESSMENT SUMMARY REPORT")
        report.append("="*80)
        
        # Building information
        report.append("\n📊 PROJECT INFORMATION")
        report.append("-"*40)
        report.append(f"Building Type: {self.building_info.get('building_type', 'N/A')}")
        report.append(f"Performance Level: {self.building_info.get('performance_level', 'N/A')}")
        report.append(f"Number of Storeys: {self.building_info.get('num_storeys', 'N/A')}")
        report.append(f"Total Height: {self.building_info.get('total_height', 'N/A')} m")
        report.append(f"Total Seismic Weight: {self.base_shear_results['total_weight']:,.0f} kN")
        
        # Seismic hazard
        report.append("\n🌋 SEISMIC HAZARD")
        report.append("-"*40)
        report.append(f"Site Class: {self.building_info.get('site_class', 'N/A')}")
        report.append(f"S_s: {self.building_info.get('S_s', 'N/A')} g")
        report.append(f"S_1: {self.building_info.get('S_1', 'N/A')} g")
        report.append(f"Seismicity Level: {self.hazard_params['seismicity_level']}")
        report.append(f"Building Period: {self.base_shear_results['T']:.3f} s")
        report.append(f"Spectral Acceleration: {self.base_shear_results['S_a']:.3f} g")
        
        # Base shear
        report.append("\n⚖️ BASE SHEAR RESULTS")
        report.append("-"*40)
        report.append(f"Base Shear, V: {self.base_shear_results['base_shear']:,.0f} kN")
        report.append(f"Base Shear Coefficient: {self.base_shear_results['S_a_T1']:.3f}")
        
        # Structural checks
        report.append("\n🏗️ STRUCTURAL CHECKS SUMMARY")
        report.append("-"*40)
        
        if self.checks:
            checks_df = pd.DataFrame(self.checks)
            report.append(checks_df.to_string(index=False))
        else:
            report.append("No structural checks performed.")
        
        # Tier 1 procedure
        report.append("\n📋 TIER 1 PROCEDURE")
        report.append("-"*40)
        tier1_allowed = self.check_tier1_limitations()
        if tier1_allowed:
            report.append("Status: ✅ Tier 1 procedure is permitted for this building")
        else:
            report.append("Status: ❌ Tier 1 procedure is NOT permitted - consider Tier 2 or 3")
        
        report.append("\n" + "="*80)
        report.append("END OF REPORT")
        report.append("="*80)
        
        return "\n".join(report)
    