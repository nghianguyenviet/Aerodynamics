# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:28:05 2023

@author: nghianguyen
"""
'''
** Definition
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
        