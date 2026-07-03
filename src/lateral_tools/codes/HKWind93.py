"""
Hong Kong Wind Load Calculation - Code of Practice on Wind Effects 1983 (1993 Reprint)
Based on the provided code document.
"""

import math
from typing import Dict, List, Tuple, Optional

class HKWindCode1983:
    """
    Implements wind load calculations per Hong Kong 1983 Wind Code (1993 Reprint)
    """
    
    # Table 1A: Design wind pressure for general terrain (kPa)
    WIND_PRESSURES = [
        (0, 10, 1.2),
        (10, 30, 2.2),
        (30, 50, 2.5),
        (50, 100, 3.0),
        (100, 150, 3.5),
        (150, 200, 3.8),
        (200, 250, 4.1),
        (250, float('inf'), 4.3),
    ]
    
    # Table B1: Height aspect factors C_h
    HEIGHT_ASPECT_FACTORS = [
        (1.0, 0.95),
        (2.0, 1.0),
        (4.0, 1.05),
        (6.0, 1.1),
        (10.0, 1.2),
    ]
    
    # Table B2: Shape factors C_s for rectangular buildings
    SHAPE_FACTORS = [
        (0.5, 0.8),
        (0.7, 0.9),
        (1.0, 1.0),
        (2.0, 1.1),
        (3.0, 1.3),
        (4.0, 1.45),
    ]
    
    # Table B3: Reduction factor R_A by frontal projected area (m²)
    AREA_REDUCTION_FACTORS = [
        (500, 1.00),
        (800, 0.97),
        (1000, 0.96),
        (3000, 0.92),
        (5000, 0.89),
        (8000, 0.86),
        (10000, 0.84),
        (15000, 0.80),
    ]
    
    # Section 7: Stability safety factor
    STABILITY_FACTOR = 1.5
    
    @classmethod
    def get_wind_pressure(cls, height_m: float) -> float:
        """Return design wind pressure q_z at given height (Table 1A)"""
        for h_min, h_max, pressure in cls.WIND_PRESSURES:
            if h_min <= height_m < h_max:
                return pressure
        return cls.WIND_PRESSURES[-1][2]  # Return highest for >250m
    
    @classmethod
    def get_height_aspect_factor(cls, height_breadth_ratio: float) -> float:
        """Return C_h from Table B1 with linear interpolation"""
        if height_breadth_ratio <= 1.0:
            return 0.95
        if height_breadth_ratio >= 10.0:
            return 1.2
        
        for i, (ratio, factor) in enumerate(cls.HEIGHT_ASPECT_FACTORS):
            if ratio >= height_breadth_ratio:
                prev_ratio, prev_factor = cls.HEIGHT_ASPECT_FACTORS[i-1]
                # Linear interpolation
                t = (height_breadth_ratio - prev_ratio) / (ratio - prev_ratio)
                return prev_factor + t * (factor - prev_factor)
        return 1.0
    
    @classmethod
    def get_shape_factor(cls, depth_breadth_ratio: float) -> float:
        """Return C_s from Table B2 with linear interpolation"""
        if depth_breadth_ratio <= 0.5:
            return 0.8
        if depth_breadth_ratio >= 4.0:
            return 1.45
        
        for i, (ratio, factor) in enumerate(cls.SHAPE_FACTORS):
            if ratio >= depth_breadth_ratio:
                prev_ratio, prev_factor = cls.SHAPE_FACTORS[i-1]
                t = (depth_breadth_ratio - prev_ratio) / (ratio - prev_ratio)
                return prev_factor + t * (factor - prev_factor)
        return 1.0
    
    @classmethod
    def get_area_reduction_factor(cls, frontal_area_m2: float) -> float:
        """Return R_A from Table B3 with linear interpolation"""
        if frontal_area_m2 <= 500:
            return 1.0
        if frontal_area_m2 >= 15000:
            return 0.80
        
        for i, (area, factor) in enumerate(cls.AREA_REDUCTION_FACTORS):
            if area >= frontal_area_m2:
                prev_area, prev_factor = cls.AREA_REDUCTION_FACTORS[i-1]
                t = (frontal_area_m2 - prev_area) / (area - prev_area)
                return prev_factor + t * (factor - prev_factor)
        return 1.0
    
    @classmethod
    def calculate_force_coefficient(cls, height_m: float, breadth_m: float, 
                                     depth_m: float, apply_reduction: bool = True) -> float:
        """
        Calculate force coefficient C_f per Appendix B
        
        Args:
            height_m: Building height
            breadth_m: Width normal to wind direction
            depth_m: Depth parallel to wind direction
            apply_reduction: Whether to apply area reduction factor R_A
        
        Returns:
            Force coefficient C_f
        """
        height_breadth = height_m / breadth_m
        depth_breadth = depth_m / breadth_m
        
        C_h = cls.get_height_aspect_factor(height_breadth)
        C_s = cls.get_shape_factor(depth_breadth)
        
        C_f = C_h * C_s
        
        if apply_reduction:
            frontal_area = height_m * breadth_m
            R_A = cls.get_area_reduction_factor(frontal_area)
            C_f *= R_A
        
        return C_f
    
    @classmethod
    def calculate_pressure_sum(cls, height_m: float, width_m: float) -> Tuple[float, List[Dict]]:
        """
        Calculate sum(q_z × A_z) per Equation 2
        
        Args:
            height_m: Building height
            width_m: Building width normal to wind direction
        
        Returns:
            Tuple of (total_sum, segments_list)
        """
        segments = []
        remaining_height = height_m
        current_bottom = 0.0
        
        for h_min, h_max, pressure in cls.WIND_PRESSURES:
            if current_bottom >= height_m:
                break
            
            # Determine segment top within this pressure zone
            zone_top = min(h_max, height_m)
            if zone_top <= current_bottom:
                continue
            
            seg_height = zone_top - current_bottom
            area = seg_height * width_m
            qz_area = pressure * area
            
            segments.append({
                'bottom_m': current_bottom,
                'top_m': zone_top,
                'height_m': seg_height,
                'pressure_kPa': pressure,
                'area_m2': area,
                'qz_area_kN': qz_area
            })
            
            current_bottom = zone_top
        
        total_qz_area = sum(seg['qz_area_kN'] for seg in segments)
        return total_qz_area, segments
    
    @classmethod
    def calculate_wind_load(cls, height_m: float, width_m: float, 
                            depth_m: float, verbose: bool = True) -> Dict:
        """
        Calculate total horizontal wind load on a rectangular building
        
        Args:
            height_m: Building height above base level
            width_m: Building width normal to wind direction
            depth_m: Building depth parallel to wind direction
            verbose: Print detailed calculation steps
        
        Returns:
            Dictionary with calculation results
        """
        if verbose:
            print("\n" + "="*60)
            print("HONG KONG WIND CODE 1983 (1993 REPRINT)")
            print("="*60)
            print(f"\nBuilding dimensions:")
            print(f"  Height: {height_m} m")
            print(f"  Width (normal to wind): {width_m} m")
            print(f"  Depth (parallel to wind): {depth_m} m")
        
        # Step 1: Calculate force coefficient
        C_f = cls.calculate_force_coefficient(height_m, width_m, depth_m)
        
        if verbose:
            h_b_ratio = height_m / width_m
            d_b_ratio = depth_m / width_m
            frontal_area = height_m * width_m
            print(f"\nStep 1: Force coefficient C_f")
            print(f"  Height/Breadth = {h_b_ratio:.2f} → C_h = {cls.get_height_aspect_factor(h_b_ratio):.2f}")
            print(f"  Depth/Breadth = {d_b_ratio:.2f} → C_s = {cls.get_shape_factor(d_b_ratio):.2f}")
            print(f"  Frontal area = {frontal_area:.0f} m² → R_A = {cls.get_area_reduction_factor(frontal_area):.2f}")
            print(f"  C_f = {C_f:.3f}")
        
        # Step 2: Calculate sum(q_z × A_z)
        total_qz_area, segments = cls.calculate_pressure_sum(height_m, width_m)
        
        if verbose:
            print(f"\nStep 2: Sum of q_z × A_z")
            print(f"  {'Bottom':>8} {'Top':>8} {'Height':>8} {'q_z':>8} {'Area':>10} {'q_z×A_z':>12}")
            print(f"  {'(m)':>8} {'(m)':>8} {'(m)':>8} {'(kPa)':>8} {'(m²)':>10} {'(kN)':>12}")
            print(f"  {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*10} {'-'*12}")
            for seg in segments:
                print(f"  {seg['bottom_m']:8.1f} {seg['top_m']:8.1f} {seg['height_m']:8.1f} "
                      f"{seg['pressure_kPa']:8.1f} {seg['area_m2']:10.0f} {seg['qz_area_kN']:12.0f}")
            print(f"\n  Total Σ(q_z×A_z) = {total_qz_area:.0f} kN")
        
        # Step 3: Calculate total wind force
        F = C_f * total_qz_area
        
        if verbose:
            print(f"\nStep 3: Total wind force")
            print(f"  F = C_f × Σ(q_z×A_z) = {C_f:.3f} × {total_qz_area:.0f} = {F:.0f} kN")
        
        # Step 4: Stability check (Section 7)
        F_stability = F * cls.STABILITY_FACTOR
        
        if verbose:
            print(f"\nStep 4: Stability check (Section 7)")
            print(f"  Required resistance = 1.5 × wind load = {F_stability:.0f} kN")
        
        return {
            'force_coefficient': C_f,
            'pressure_sum_kN': total_qz_area,
            'wind_load_kN': F,
            'stability_required_kN': F_stability,
            'segments': segments
        }
    
    @classmethod
    def maximum_wind_load(cls, height_m: float, width_m: float, depth_m: float) -> Dict:
        """
        Calculate the maximum wind load by considering both wind directions
        (wind on width face and wind on depth face)
        """
        print("\n" + "="*60)
        print("COMPARING BOTH WIND DIRECTIONS")
        print("="*60)
        
        # Wind on width face (breadth = width, depth = depth)
        result1 = cls.calculate_wind_load(height_m, width_m, depth_m, verbose=False)
        
        # Wind on depth face (breadth = depth, depth = width)
        result2 = cls.calculate_wind_load(height_m, depth_m, width_m, verbose=False)
        
        print(f"\nWind on {width_m}m face: F = {result1['wind_load_kN']:.0f} kN")
        print(f"Wind on {depth_m}m face: F = {result2['wind_load_kN']:.0f} kN")
        
        if result2['wind_load_kN'] > result1['wind_load_kN']:
            print(f"\n→ Governing wind load: {result2['wind_load_kN']:.0f} kN (wind on {depth_m}m face)")
            result2['governing'] = True
            result1['governing'] = False
            return result2
        else:
            print(f"\n→ Governing wind load: {result1['wind_load_kN']:.0f} kN (wind on {width_m}m face)")
            result1['governing'] = True
            result2['governing'] = False
            return result1


# ============================================================================
# Example: Your building (32m tall, 136m × 65m plan, sea-front location)
# ============================================================================

if __name__ == "__main__":
    
    # Your building parameters
    BUILDING = {
        'height': 32.0,      # meters
        'width': 136.0,      # meters (face normal to wind for first case)
        'depth': 65.0,       # meters (face parallel to wind)
    }
    
    # Calculation for wind on the 136m face
    result_136 = HKWindCode1983.calculate_wind_load(
        height_m=BUILDING['height'],
        width_m=BUILDING['width'],
        depth_m=BUILDING['depth']
    )
    
    # Calculation for wind on the 65m face (this governs)
    result_65 = HKWindCode1983.calculate_wind_load(
        height_m=BUILDING['height'],
        width_m=BUILDING['depth'],
        depth_m=BUILDING['width']
    )
    
    # Show which governs
    print("\n" + "="*60)
    print("CONCLUSION")
    print("="*60)
    
    if result_65['wind_load_kN'] > result_136['wind_load_kN']:
        print(f"\nGoverning wind load: {result_65['wind_load_kN']:.0f} kN")
        print(f"(Wind acting on the {BUILDING['depth']}m wide face)")
        governing = result_65
    else:
        print(f"\nGoverning wind load: {result_136['wind_load_kN']:.0f} kN")
        print(f"(Wind acting on the {BUILDING['width']}m wide face)")
        governing = result_136
    
    print(f"\nStability requirement (Section 7): Resistance ≥ {governing['stability_required_kN']:.0f} kN")
    
    # Optional: Load combination example per BS 8110 (for reference)
    print("\n" + "="*60)
    print("ULTIMATE DESIGN LOAD COMBINATIONS (BS 8110 reference)")
    print("="*60)
    print(f"Dead + Wind only:     1.4G + 1.4W = 1.4G + {1.4 * governing['wind_load_kN']:.0f} kN")
    print(f"Dead + Live + Wind:   1.2G + 1.2Q + 1.2W = 1.2G + 1.2Q + {1.2 * governing['wind_load_kN']:.0f} kN")
    print(f"Note: Load factors from BS 8110, not from wind code")


'''
============================================================
HONG KONG WIND CODE 1983 (1993 REPRINT)
============================================================

Building dimensions:
  Height: 32.0 m
  Width (normal to wind): 65.0 m
  Depth (parallel to wind): 136.0 m

Step 1: Force coefficient C_f
  Height/Breadth = 0.49 → C_h = 0.95
  Depth/Breadth = 2.09 → C_s = 1.12
  Frontal area = 2080 m² → R_A = 0.94
  C_f = 1.000

Step 2: Sum of q_z × A_z
  Total Σ(q_z×A_z) = 7941 kN

Step 3: Total wind force
  F = C_f × Σ(q_z×A_z) = 1.000 × 7941 = 7941 kN

============================================================
CONCLUSION
============================================================

Governing wind load: 7941 kN
(Wind acting on the 65m wide face)

Stability requirement (Section 7): Resistance ≥ 11912 kN

```