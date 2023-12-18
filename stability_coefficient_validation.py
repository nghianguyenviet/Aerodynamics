# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:28:05 2023

@author: nghianguyen
"""
'''
** Definition

----------------- Longditudinal & Lateral Stability --------------------

CD_0: aircraft drag coef. at equilibrium flight condition 
CL_0: aircraft lift coef. at equilibrium flight condition 
CL_alpha: (other denote: a) lift curve slope of aircraft ≈ lift curve slope of wing, rad
theta_0 : angle of rotation, at equilibrium flight condition, of X_A axis from horizontal plane, rad or deg
wing_aspect_ratio:  Aspect ratio of wing
aircraft_effi_factor : airplane efficiency factor
effi_factor_HT : horizontal tail efficiency factor
l_T : distance, wing quarter chord to horizontal tail quarter chord, ft
area_HT : horizontal tail area (including elevator), ft^2:
S : gross wing area, ft^2
CL_alpha_HT : (other denote: CL_alphCL_alpha_HT) ,lift curve slope of horizontal tail, rad^-1 
downwash_rate_alpha : (denote: dϵ/dα) rate of change of downwash angle at tail with change in angle of attack, rad/rad
c : mean aerodynamic chord
d_alphCL_alpha_HT_by_delta_e : elevator effectiveness
stick_fixed_static_margin : (denote: h - h_n) stick-fixed static margin, fraction of c

area_VT : vertical tail area (including rudder), ft^2
CL_alpha_VT : lift curve slope of vertical tail, deg^-1 or rad^-1
downwash_rate_beta : rate of change of sidewash angle at vertical tail with change in sideslip angle, rad/rad
span_VT : vertical tail span, ft
lambda_F : tape ratio of vertical tail
l_F : distance, wing quarter chord to vertical tail quarter chord, ft
rudder_effectiveness : rudder effectiveness
vol_fuse : volume of fuselage, ft^3
h : fuselage height at wing root
w : fuselage width at wing root
b : wing span
CD_w : wing drag coef.
rudder_effectiveness :rudder effectiveness
Cl_beta_w_by_r1 : page 19
Cl_beta_w_by_r2 : page 19
aileron_effectiveness : aileron effectiveness
Cl_p : page 20
Z_w : vertical distance of root quarter chord below fuselage center line, ft
Z_F : vertival distance, aircraft CG to mean chord of vertival tail, ft
y1 : distance, aircraft centerline to inboard end of ailereon, ft
y2 : distance, aircraft centerline to outboard end of aileron, ft

# Note: All of the equations are taken from LINEARIZED MATHEMATICAL MODELS FOR DE HAVILLAND CANADA "BUFFALO & TWIN OTTER" STOL TRANSPORTS
'''
# VT: Vertical Tail
# HT: Horizontal Tail
# l_T : distance, wing quarter chord to horizontal tail quarter chord, ft
# l_F : distance, wing quarter chord to vertical tail quarter chord, ft
# Z_w : vertical distance of root quarter chord below fuselage center line, ft
# Z_F : vertival distance, aircraft CG to mean chord of vertival tail, ft
# Cl_beta_w_by_r1 : page 19 in the ref
# Cl_beta_w_by_r2 : page 19 in the ref
# y1 : distance, aircraft centerline to inboard end of ailereon, ft
# y2 : distance, aircraft centerline to outboard end of aileron, ft
# theta_0: angle of rotation, at equilibrium flight condition, of X_A axis from horizontal plane, rad or deg

###### Dynamic Stability Coefficients Estimating ########################
import numpy as np

def calc_cz_q(CL_alpha_HT,effi_factor_HT, l_T, area_HT, mac, area_main_wing):
    return -2 * CL_alpha_HT * horizontal_tail_vol(effi_factor_HT,l_T, area_HT, mac, area_main_wing) 

def calc_cm_q(l_T, mac,CL_alpha_HT, effi_factor_HT, area_HT, area_main_wing):
    return calc_cz_q(CL_alpha_HT, effi_factor_HT, l_T, area_HT, mac, area_main_wing) * l_T / mac 

def calc_cy_p(CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT):
    return -2/3 * CL_alpha_VT * (span_VT/wing_span) * (area_VT/area_main_wing) * ((1+ 2*tape_ratio_VT)/(1+tape_ratio_VT)) 


def calc_cn_p (CL_0,CL_alpha, wing_aspect_ratio,CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT,l_F):
    return (-CL_0/4 * (1 - CL_alpha/(np.pi*wing_aspect_ratio)) 
            -calc_cy_p(CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT) * l_F /wing_span)

def calc_cn_r(l_F, wing_span, CD_w, CL_alpha_VT, area_VT, area_main_wing):
    return -(l_F/wing_span) * calc_cy_r(CL_alpha_VT, l_F, area_VT, wing_span, area_main_wing) - CD_w/4 

def calc_cl_p (cl_p):
    return cl_p 
    
def calc_cl_r (CL_0, Z_F, wing_span,CL_alpha_VT,l_F,area_VT,area_main_wing):
    return CL_0 / 4 + calc_cy_r(CL_alpha_VT, l_F, area_VT, wing_span, area_main_wing) * Z_F/wing_span 

def horizontal_tail_vol(effi_factor_HT, l_T , area_HT, mac, area_main_wing):
    return effi_factor_HT * (l_T / mac) * (area_HT / area_main_wing)

def calc_cz_alpha_dot(CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, l_T, area_HT, mac, area_main_wing):
    return -2 * CL_alpha_HT * horizontal_tail_vol(effi_factor_HT,l_T, area_HT, mac, area_main_wing) * downwash_rate_alpha

def calc_cm_alpha_dot(l_T, mac, CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, area_HT, area_main_wing): 
    return calc_cz_alpha_dot(CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, l_T, area_HT, mac, area_main_wing) * l_T / mac

######################################### Static Stability Coefficients Estimating #######################################################

def calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT):
    return CL_alpha_HT * elevator_effectiveness * (area_HT / area_main_wing) *effi_factor_HT

def calc_cm_delta_e(l_T, mac,CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT):
    return calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT) * l_T / mac

def calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness):
    return -CL_alpha_VT * (area_VT/area_main_wing) * rudder_effectiveness

def calc_cn_sigma_a():
    return 0

def calc_cn_sigma_r(l_F,wing_span,CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness):
    return -calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness) * l_F/wing_span

def calc_cl_beta (Cl_beta_w_by_r1, wing_dihedral_angle, Cl_beta_w_by_r2,area_VT, area_main_wing, 
             CL_alpha_VT, downwash_rate_beta,wing_aspect_ratio,Z_w,wing_span,fuselage_height,fuselage_width,Z_F):
    
    return ((Cl_beta_w_by_r1 - Cl_beta_w_by_r2) * wing_dihedral_angle + 1.2*np.sqrt(wing_aspect_ratio) * (Z_w/wing_span) 
            * (fuselage_height + fuselage_width)/wing_span + calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta)*Z_F/wing_span)

def calc_cl_sigma_a(CL_alpha, aileron_effectiveness, area_main_wing, wing_span, y1, y2):
    return (2 * CL_alpha * aileron_effectiveness / (area_main_wing * wing_span)) * (7.41*(y2**2 - y1**2) - 0.062 * (47**3 - 33**3))

def calc_cl_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness, wing_span, Z_F):
    return calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness) * Z_F/wing_span


def calc_cx_u(CD_0,CL_0,theta_0 = 0):
    return - 3 * CD_0 - CL_0 * np.tan(theta_0) 

def calc_cx_alpha(CL_0, CL_alpha, wing_aspect_ratio, aircraft_effi_factor):
    return CL_0 - 2* CL_0 * CL_alpha / (np.pi * aircraft_effi_factor * wing_aspect_ratio)

def calc_cz_u(CL_0):
    return -2 * CL_0

def calc_cz_alpha(CL_alpha, CD_0):
    return -(CL_alpha + CD_0)

def calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta):
    return -CL_alpha_VT * area_VT / area_main_wing * (1 - downwash_rate_beta)

def calc_cy_r(CL_alpha_VT, l_F, area_VT, wing_span, area_main_wing):
        return 2 * CL_alpha_VT * (l_F/wing_span) * (area_VT/area_main_wing)

def calc_cn_beta(l_F, wing_span, vol_fuse,area_main_wing,fuselage_height,fuselage_width,area_VT, CL_alpha_VT, downwash_rate_beta):
    return (-calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta) * l_F/wing_span 
            - 1.3 * vol_fuse/(area_main_wing*wing_span) * fuselage_height/fuselage_width)




effi_factor_HT = 1.0
l_T = 46
area_HT = 233
mac = 10.1
area_main_wing = 945
CL_alpha_HT = 3.5
CL_alpha_VT = 2.5
span_VT = 14
wing_span = 96
area_VT =152
tape_ratio_VT = 0.64
l_F = 44
CD_w = 0.025
CL_0 = 0.77
CD_0 = 0.057
wing_aspect_ratio = 9.75
Z_F = 10
Z_w = -3.5
y1 = 33
y2 = 47
cl_p = -0.53
CL_alpha = 5.2
downwash_rate_alpha = 0.17
elevator_effectiveness = 0.54
aileron_effectiveness = 0.43
rudder_effectiveness = 0.58
theta_0 = 0 
downwash_rate_beta = 0.1
vol_fuse = 4600
fuselage_height = 9.8
fuselage_width = 10
aircraft_effi_factor = 0.75
wing_dihedral_angle = 0.087
Cl_beta_w_by_r1 = -0.84
Cl_beta_w_by_r2 = -0.14



cz_q = calc_cz_q(CL_alpha_HT, effi_factor_HT, l_T, area_HT, mac, area_main_wing)
cm_q = calc_cm_q(l_T, mac,CL_alpha_HT, effi_factor_HT, area_HT, area_main_wing)
cy_p = calc_cy_p(CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT)
cy_r = calc_cy_r(CL_alpha_VT, l_F, area_VT, wing_span, area_main_wing)
cn_p = calc_cn_p(CL_0,CL_alpha, wing_aspect_ratio,CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT,l_F)
cn_r = calc_cn_r(l_F, wing_span, CD_w, CL_alpha_VT, area_VT, area_main_wing)
cl_p_value = calc_cl_p(cl_p)
cl_r = calc_cl_r(CL_0, Z_F, wing_span,CL_alpha_VT,l_F,area_VT,area_main_wing)
cz_alpha_dot = calc_cz_alpha_dot(CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, l_T, area_HT, mac, area_main_wing)
cz_delta_e = calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT)
cm_alpha_dot = calc_cm_alpha_dot(l_T, mac, CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, area_HT, area_main_wing)
cm_delta_e = calc_cm_delta_e(l_T, mac,CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT)
cy_sigma_r = calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness)
cn_sigma_a = calc_cn_sigma_a()
cn_sigma_r = calc_cn_sigma_r(l_F,wing_span,CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness)
cl_sigma_a = calc_cl_sigma_a(CL_alpha, aileron_effectiveness, area_main_wing, wing_span, y1, y2)
cl_sigma_r = calc_cl_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness, wing_span, Z_F)
cx_u = calc_cx_u(CD_0, CL_0, theta_0)
cx_alpha = calc_cx_alpha(CL_0,CL_alpha, wing_aspect_ratio, aircraft_effi_factor)
cz_u = calc_cz_u(CL_0)
cz_alpha = calc_cz_alpha(CL_alpha,CD_0)
cy_beta = calc_cy_beta(area_VT,area_main_wing,CL_alpha_VT,downwash_rate_beta)
cn_beta = calc_cn_beta(l_F,wing_span,vol_fuse,area_main_wing,fuselage_height,fuselage_width,area_VT,CL_alpha_VT,downwash_rate_beta) 
cl_beta = calc_cl_beta (Cl_beta_w_by_r1, wing_dihedral_angle, Cl_beta_w_by_r2,area_VT, area_main_wing, 
             CL_alpha_VT, downwash_rate_beta,wing_aspect_ratio,Z_w,wing_span,fuselage_height,fuselage_width,Z_F)

print("------------------------Calculation of Dynamic stability coefficients -----------------------\n")
print("cz_q:", round(cz_q, 3))
print("cm_q:", round(cm_q, 3))
print("cy_p:", round(cy_p, 3))
print("cy_r:", round(cy_r, 3))
print("cn_p:", round(cn_p, 3))
print("cn_r:", round(cn_r, 3))
print("cl_p:", round(cl_p_value, 3))
print("cl_r:", round(cl_r, 3))
print("cz_alpha_dot:", round(cz_alpha_dot, 3))
print("cm_alpha_dot:", round(cm_alpha_dot, 3))
print('\n')

# print("------------------------Calculation of Static stability coefficients -----------------------\n")

# print("cz_delta_e:", round(cz_delta_e, 3))
# print("cm_delta_e:", round(cm_delta_e, 3))
# print("cy_sigma_r:", round(cy_sigma_r, 3))
# print("cn_sigma_a:", round(cn_sigma_a, 3))
# print("cn_sigma_r:", round(cn_sigma_r, 3))
# print("cl_sigma_a:", round(cl_sigma_a, 3))
# print("cl_sigma_r:", round(cl_sigma_r, 3))
# print("cx_u =", round(cx_u, 3))
# print("cx_alpha =", round(cx_alpha, 3))
# print("cz_u =", round(cz_u, 3))
# print("cz_alpha =", round(cz_alpha, 3))
# print("cy_beta =", round(cy_beta, 3))
# print("cn_beta =", round(cn_beta, 3))
# print("cl_beta =", round(cl_beta, 3))
