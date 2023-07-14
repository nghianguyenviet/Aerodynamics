# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:28:05 2023

@author: nghianguyen
"""
'''
** Definition

----------------- Longditudinal & Lateral Stability --------------------

AR:  Aspect ratio
CD_0: Reference drag coefficient
CL_0: Reference lift coefficient
CL_alpha: Airplane lift curve slope
CL_alpha_w: Wing lift curve slope
CL_alpha_t: Tail lift curve slope
c_mean: Mean aerodynamic chord
e_factor: Osward's span efficiency factor
l_t: Distance from center of gravity to tail quarter chord
V_H: Horizontal tail volume ratio
M: Flight mach number
S: Wing area
S_t: Horizontal tail area
d_epsi_by_d_alpha: Change in downwash due to a change in angle of attack
eta: Efficiency factor of the horizontal tail
d_CL_t: ??????????
d_delta_e: ??????????
M_0: ????
Cm_alpha_fus: ?????
c: ????
b: Wing span
K: Empirical factor
l_v: Distance from center of gravity to vertical tail aerodynamic center
V_v: Vertical tail volume ratio
S_v: Vertival tail area
z_v: Distance from center of pressure of vertical tail to fuselage centerline
wing_dihedral_angle
wing_sweep_angle
eta_v = Efficiency factor of the vertical tail
taper_ratio (tip chord/root chord)
d_sigma_by_d_beta: Change in sidewash angle with a change in sideslip angle
tau : ?????????
Cy_beta_tail: ?????????
integral_cy_dy: ????????

'''

import numpy as np

class longitudinal_stability_coefficients:
    
    # X-force derivatives
    
    def calc_CX_u(CD_u, CD_0, CT_u):
        return -(CD_u + 2 * CD_0) + CT_u
    
    def calc_CX_alpha(CL_0, Cl_alpha, AR, e_factor):
        return CL_0 - 2* CL_0 * Cl_alpha / (np.pi * e_factor * AR)
    
    # Z-force derivatives
    
    def calc_Cz_u(M, CL_0):
        return -M**2 * CL_0 / (1 - M**2) + 2*CL_0
    
    def calc_CZ_alpha(CL_alpha, CD_0):
        return -(CL_alpha + CD_0)
    
    def calc_CZ_alpha_dot(eta, CL_alpha_t, V_H, d_epsi_by_d_alpha):
        return -2 * eta * CL_alpha_t * V_H * d_epsi_by_d_alpha
    
    def calc_Cz_q(eta, CL_alpha_t, V_H):
        return -2 * eta * CL_alpha_t * V_H
    
    def calc_CZ_delta_e(S_t, S, eta, d_CL_t, d_delta_e):
        return -S_t * eta *d_CL_t / (S * d_delta_e)
    
    def calc_CL_delta_e(S_t, S, eta, d_CL_t_by_d_delta_e):
        return S_t * eta * d_CL_t_by_d_delta_e / S
    
    # Pitching moment derivatives
    
    def calc_Cm_u(d_Cm_by_d_M, M0):
        return d_Cm_by_d_M * M0
    
    def calc_Cm_alpha(CL_alpha_w, X_cg, X_ac, c_mean, Cm_alpha_fus, eta, V_H, CL_alpha_t, d_epsi_by_d_alpha):
        return CL_alpha_w * (X_cg / c_mean - X_ac / c_mean) + Cm_alpha_fus - eta * V_H * CL_alpha_t * (1 - d_epsi_by_d_alpha)
    
    def calc_Cm_alpha_dot(eta, CL_alpha_t, V_H, l_t, c, d_epsi_by_d_alpha):
        return -2 * eta * CL_alpha_t * V_H * l_t * d_epsi_by_d_alpha / c
    
    def calc_Cm_q(eta, CL_alpha_t, V_H, l_t, c):
        return -2 * eta * CL_alpha_t * V_H * l_t / c
    
    def calc_Cm_delta_e(eta, V_H, d_Cl_t_by_d_delta_e):
        return -eta * V_H * d_Cl_t_by_d_delta_e
    
class lateral_stability_coefficients:
    
    # Y-force derivatives
    
    def Cy_beta(eta, S_v, S, CL_alpha_v, d_sigma_by_d_beta):
        return -eta * S_v * CL_alpha_v * (1 + d_sigma_by_d_beta) / S
    
    def Cy_p(CL, AR, wing_sweep_angle):
        return CL * np.tan(wing_sweep_angle) * (AR + np.cos(wing_sweep_angle)) / (AR + 4 * np.cos(wing_sweep_angle))
    
    def Cy_r(l_v, b, Cy_beta_tail):
        return -2 * l_v * Cy_beta_tail / b
    
    def Cy_sigma_r(S_v, S, tau, CL_alpha_v):
        return S_v * tau * CL_alpha_v / S
    
    # Yawing moment derivatives
    
    def Cn_beta(Cn_beta_wf, eta_v, V_v, CL_alpha_v, d_sigma_by_d_beta):
        return Cn_beta_wf + eta_v * V_v * CL_alpha_v * (1 + d_sigma_by_d_beta)
    
    def Cn_p (CL):
        return -CL / 8
    
    def Cn_r(eta_v, V_v, l_v, CL_alpha_v, b):
        return -2 * eta_v * V_v * l_v * CL_alpha_v / b
    
    def Cn_sigma_a(K, CL_0, Cl_sigma_a):
        return 2 * K * CL_0 * Cl_sigma_a
    
    def Cn_sigma_r(V_v, eta_v, tau, CL_alpha_v):
        return -V_v * eta_v * tau * CL_alpha_v
    
    # Rolling monment derivatives
    
    def Cl_beta (Cl_beta_by_wing_dihedral_angle, wing_dihedral_angle, delta_Cl_beta):
        return Cl_beta_by_wing_dihedral_angle * wing_dihedral_angle + delta_Cl_beta
    
    def Cl_p (CL_alpha, taper_ratio):
        return - CL_alpha * (1 + 3 * taper_ratio) / (12 * (1 + taper_ratio))
        
    def Cl_r (CL, l_v, z_v, Cy_beta_tail, b):
        return CL / 4 + -2 * (l_v / b) * (z_v / b) * Cy_beta_tail
    
    def Cl_sigma_a(CL_alpha, tau, integral_cy_dy, S, b):
        return 2 * CL_alpha * tau * integral_cy_dy / (S * b)
    
    def Cl_sigma_r(S_v, z_v, tau, CL_alpha_w, S, b):
        return S_v * (z_v / b) * tau * CL_alpha_w / S
        
    
    
    
    
        