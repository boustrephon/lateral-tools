""" 
Macau Regulamento de Segurança e Acções em 
Estruturas de Edifícios e Pontes  (RSAEEP)
 
Seismic Action Revision (Consultation Document)
Version - 15th of April, 2008

There are two explicit performance levels:

* No-collapse requirement
* Damage limitation requirement

There are four categories of structures based on socio-economic importance, which have associated importance factors

1. High priority - civil protection etc
2. High consequence - public assembly etc
3. Ordinary - not in other categories
4. Minor - minor importance for public safety and temporary structures

Two seismic actions are defined:

* a design seismic action (for collapse prevention)  - 10% probability of exceedance in 50 years (475 years)
* a serviceability seismic action (for damage limitation) - 10% probability of exceedance in 10 years (95 years)

The $\nu$ reduction factor is a multiplier on the reference (characteristic) seismic action to generate the serviceability seismic action, equivalent to decreasing the return factor from 475 to 95 years.

- *Building* - construction facility basically to be used directly by people; 
- *Return period or average recurrence interval of an event* - time period (in 
years) equal to the inverse of its annual probability of exceedance (expressed 
in decimal format), or equal to its probability of exceedance (in %) during 
the life-time (or reference period) of a structure. For example: an event with 
a return period of 475 years, will occur in average once in every 475 years, 
and will have an annual probability to be exceeded of 0.002, or 10% of 
probability of exceedance during an intended life-time of a structure of 50 
years; 
- *China Seismic Ground Zonation Maps 1 and 2* - official document where it is 
assigned, also for the Macau region, the seismic input for the following 
seismic hazard: peak ground acceleration (PGA) value and the spectrum 
characteristic period (Tg), for a ground site-class type II, with a 10% 
probability of exceedance in 50 years (equivalent to a return period of 475 
years); 
- *Reference or characteristic peak ground acceleration (PGA)* - defines the 
seismic hazard on type II site-class, and it is derived, for Macau region, from 
the “China Zonation Map” but adjusted for local geological conditions. It 
corresponds to the characteristic return period of 475 years, and it is the 
basis for the determination of the design seismic actions, to be considered in 
the analysis of structures qualified in ordinary importance category;  
- Spectrum characteristic period, Tg - spectral period parameter of the response 
spectrum curve, depending of the ground site-class type, where the spectrum 
curve starts its decay section; 
- *Importance category of a building or structure* - qualification of structures 
according to the social-economic importance of their using functions. Higher 
importance levels should correspond to lower severity grades of after-shake 
damaging effects. Reflected through the corresponding importance factor; 
- *Acceleration response spectrum* - is a plot of the maxima of the acceleration 
response for a single degree of freedom (SDOF) system with various 
fundamental vibration periods, when subjected to an earthquake ground 
motion; 
- *Reduced-elastic or design response spectrum for non-linear-elastic response* - 
is an elastic spectrum corresponding to a target of inelasticity and is 
estimated by dividing the ordinates of the elastic spectra by the behavior 
factor (q). Is to be used in the assessment of the seismic effects in the 
ultimate (or no-local-collapse) limit state, and so called as design spectrum; 
- *Behavior factor (q)* - is the quantification of the inelastic behavior of 
structures, closely related to their ductility performance. The value q is the 
ratio of the seismic forces that a structure would experience if its response 
was full elastic with 5% damping, to the (smaller) seismic yield forces that 
could occur in the system, and so enabling some grade of inelastic 
deformation occurrence. This is equivalent to reduce the seismic elastic 
values induced in structures (or the elastic response spectrum ordinates), 
dividing them by a behavior factor (q) and then performing an explicit 
elastic analysis. It should be pointed out that high values of the behavior 
factor (q) correspond to large inelastic deformations, normally equal to the 
ones derived from the linear analysis, multiplied by the above defined 
behavior factor (q) Otherwise, if non-linear analysis methods were used, 
then both the effective internal forces and displacements are taken to be 
equal to those directly derived from the non-linear analysis, therefore 
without the intervention of the behavior factor(q); 
- *Reduced-elastic response spectrum for serviceable verification* - is the design 
response spectrum multiplied by the reduction factor (n) to be used in the 
assessment of the seismic effects in the serviceable limit state (or damage 
limitation). 
- *Reduction factor (n)* - factor to be applied to the characteristic seismic action 
in order to obtain the serviceability seismic action (for damage limitation). 

The behaviour factor (q) reflects the level of dissipation that is required
in the structure (directly correlated with ductility required) and this should
be reflected in the detailing. If `q` is equal to or greater than 3, then there is
the requirement for high dissipative structural behaviour, and the steel and
RC codes should be consulted.

* `q` l`ess than or equal to 1.5 only requires low dissipative structural behaviour (low ductility). 
* `q` less than 3 requires medium dissipative structural behaviour (medium ductility).
* `q` equal to or higher than 3 requires high dissipative structural behaviour (high ductility).

Regularity

Both plan and elevation irreularities may be present. If structures are irregular in 
elevation, then the `q`-factor shall be multiplied by 0.8
"""

def get_T_g(site_class:str, is_saturated_mud:bool=True)->float:
    """Spectrum characteristic period, Tg, based on site-class 
    from Table IV.5

    Args:
        site_class (_type_): 'I', 'II', 'III', or 'IV', based on 
            shear-wave velocity (see Table IV.4)
        is_saturated_mud (bool, optional): plasticity index > 19. 
            Defaults to False.

    Returns:
        float: Spectrum characteristic period, Tg (secs)
    """
    if site_class == 'I':
        return 0.35
    elif site_class == 'II':
        return 0.65 if is_saturated_mud else 0.45
    elif site_class == 'III':
        return 0.85 if is_saturated_mud else 0.65
    elif site_class == 'IV':
        return 1.10
    else:
        raise ValueError(f'Invalid site class: {site_class}')


def get_gamma_I(importance_category:str)->float:
    """Importance factor (Table IV.1) """
    if importance_category == 'D':
        
        return 0.8
    elif importance_category == 'C':
        return 1.0
    elif importance_category == 'B':
        return 1.2
    elif importance_category == 'A':
        return 1.4
    else:
        raise ValueError(f'Invalid importance category: {importance_category}')


def get_nu(importance_category:str)->float:
    """Reduction factor (Table IV.1)
    """
    if importance_category == 'D':
        return 0.4
    elif importance_category == 'C':
        return 0.4
    elif importance_category == 'B':
        return 0.5
    elif importance_category == 'A':
        return 0.5
    else:
        raise ValueError(f'Invalid importance category: {importance_category}')



def s_e_func(T_g:float, zeta=0.05, a_gR=0.123, gamma_I=1.0, vertical=False, default=None):
    """Return elastic response spectrum as a function of building period (T). Returns
    horizontal spectrum by default unless "vertical" is set to "True".

    Args:
        T_g (float): spectrum characteristic ground period (from Table IV.5)
        zeta (float, optional): critical damping ratio. Defaults to 0.05 (5%).
        a_gR (float, optional): Ground reference acceleration. Defaults to 0.123.
        gamma_I (float, optional): Importance factor (from Table IV.1). Defaults to 1.0.
        vertical (bool, optional): return vertical spectrum instead of lateral (horizontal). Defaults to False.
        default (_type_, optional): value to return if period exceeds 6 seconds. Defaults to None.

    Raises:
        ValueError: If period exceeds 6 seconds and default is None, otherwise
            it returns the default value if it is a float greater than zero, or else
            extends the function (with a minimum of zero).

    Returns:
        _type_: elastic response spectrum function as a function of building period (T)
    """
    beta_max = 2.4   # PGA scaling factor for plateau
    # g = 9.81 # m/s2
    gamma = 0.9 + (0.05 - zeta) / (0.5 + 5 * zeta)
    eta1 = max(0.0, 0.02 + 0.125 * (0.05 - zeta))
    eta2 = max(0.55, 1.0 + (0.05 - zeta) / (0.06 + 1.7 * zeta))
    a_g = a_gR * gamma_I
    alpha_max = a_g * beta_max # based on g
    if vertical:
        alpha_max = 0.5 * alpha_max
        T_g = 0.5 * T_g
    
    def func(T):
        s_0 = 0.42
        if T <= 0.1:
            return alpha_max * (s_0 + T / 0.1 * (eta2 - s_0))
        elif T <= T_g:
            return alpha_max * eta2
        elif T <= 5 * T_g:
            return alpha_max * eta2 * (T_g / T)**gamma
        elif T <= 6:
            return alpha_max * (eta2 * 0.2**gamma - eta1 * (T - 5 * T_g))
        else:
            if default is None:
                raise ValueError('T must be less than 6 secs')
            elif default < 0:
                return max(0, alpha_max * (eta2 * 0.2**gamma - eta1 * (T - 5 * T_g)))
            else:
                return default
    
    return func


def s_d_func(T_g, q=1.5, a_gR=0.123, gamma_I=1.0, zeta=0.05, g=1.0, vertical=False, default=None):
    """Return design response spectrum as a function of building period (T). Returns
    horizontal spectrum by default unless "vertical" is set to "True".

    A value other than 0.05 should be given for zeta (damping) only if q is set to 1.0 
        (i.e. elastic response). Behaviour factors (q) specified in the code already incorporate
        the effects of damping values other than 0.05.
     
    Args:
        T_g (float): spectrum characteristic ground period (from Table IV.5)
        q (float, optional): behaviour factor. Defaults to 1.5
        a_gR (float, optional): Ground reference acceleration. Defaults to 0.123.
        gamma_I (float, optional): Importance factor (from Table IV.1). Defaults to 1.0.
        zeta (float, optional): critical damping ratio. Defaults to 0.05 (5%).
        g (float, optional): the value of the acceleration due to gravity in the units
            of the spectral acceleration. Defaults to 1.0, i.e. units of gravity (g).
            For output in m/s2, set g to 9.81 m/s2.
        vertical (bool, optional): return vertical spectrum instead of lateral (horizontal). Defaults to False.
        default (_type_, optional): value to return if period exceeds 6 seconds. Defaults to None.

    Raises:
        ValueError: If period exceeds 6 seconds and default is None, otherwise
            it returns the default value if it is a float greater than zero, or else
            extends the function (with a minimum of zero).

    Returns:
        _type_: elastic response spectrum function as a function of building period (T)"""
    beta_max = 2.4   # PGA scaling factor for plateau
    gamma = 0.9 + (0.05 - zeta) / (0.5 + 5 * zeta)  # 0.9 for zeta = 0.05
    eta1 = max(0.0, 0.02 + 0.125 * (0.05 - zeta))   # 0.02 for zeta = 0.05
    eta2 = max(0.55, 1.0 + (0.05 - zeta) / (0.06 + 1.7 * zeta))  # 1.0 for zeta = 0.05
    a_g = a_gR * gamma_I
    alpha_max = a_g * beta_max / g
    if vertical:
        alpha_max = 0.5 * alpha_max
        T_g = 0.5 * T_g
    
    def func(T):
        T_0 = 0.1
        T_1 = T_g
        T_2 = 5 * T_g
        T_3 = 6.0
        s_0 = 0.28
        if T <= 0.1:
            return alpha_max * (s_0 + T / T_0 * (eta2 / q - s_0))
        elif T <= T_1:
            return alpha_max * eta2 / q
        elif T <= T_2:
            return alpha_max * eta2 / q * (T_1 / T)**gamma
        elif T <= T_3:
            return alpha_max / q * (eta2 * 0.2**gamma - eta1 * (T - T_2))
        else:
            if default is None:
                raise ValueError(f'T must be less than {T_3} secs')
            elif default < 0:
                return max(0, alpha_max * (eta2 * 0.2**gamma - eta1 * (T - T_2)))
            else:
                return default
    
    return func

    
