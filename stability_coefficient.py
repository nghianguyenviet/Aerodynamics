# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:28:05 2023

@author: nghianguyen & Vinh PHam
"""
'''
** Definition

----------------- Aerodynamic Variables --------------------

CD_0: aircraft drag coef. at equilibrium flight condition 

CL_0: aircraft lift coef. at equilibrium flight condition 

effi_factor_HT : Efficiency factor of horizontal tail (Qt/Qw) (0.8~1.2)

aircraft_effi_factor : airplane efficiency factor

CL_alpha_wing: (other denote: a) lift curve slope of aircraft ≈ lift curve slope of wing, rad

theta_0 : angle of rotation, at equilibrium flight condition, of X_A axis from horizontal plane, rad or deg

CL_alpha_HT :  lift curve slope of horizontal tail, rad^-1 


d_CL_alpha_HT_by_delta_e : elevator effectiveness

Cl_beta_w_by_r1 : page 19
Cl_beta_w_by_r2 : page 19

Cl_p : page 20
stick_fixed_static_margin : (denote: h - h_n) stick-fixed static margin, fraction of c
CL_alpha_VT : lift curve slope of vertical tail, deg^-1 or rad^-1

CD_w : wing drag coef.
rudder_effectiveness :rudder effectiveness

aileron_effectiveness : aileron effectiveness

downwash_rate_beta : rate of change of sidewash angle at vertical tail with change in sideslip angle, rad/rad
----------------- Geometric Variables --------------------


aspect_ratio_wing:  Aspect ratio of wing

tape_ratio_wing : tip_chord/root_chord

distance_HT : distance, wing quarter chord to horizontal tail quarter chord, ft

area_HT : horizontal tail area (including elevator), ft^2:

S : gross wing area, ft^2
mac_wing : mean aerodynamic chord of wing
mac_HT : mean aerodynamic chord of HT tails
mac_VT : mean aerodynamic chord of VT tails
area_VT : vertical tail area (including rudder), ft^2

span_VT : vertical tail span, ft
taper_ratio_VT : taper ratio of vertical tail
distance_VT : distance, wing quarter chord to vertical tail quarter chord, ft
rudder_effectiveness : rudder effectiveness
vol_fuse : volume of fuselage, ft^3
fuselage_height: fuselage height at wing root
fuselage_width: fuselage width at wing root
wing_span: wing span
Z_w : vertical distance of root quarter chord below fuselage center line, ft
Z_VT : vertival distance, aircraft CG to mean chord of vertival tail, ft
y1 : distance, aircraft centerline to inboard end of ailereon, ft
y2 : distance, aircraft centerline to outboard end of aileron, ft
# Note: All of the equations are taken from LINEARIZED MATHEMATICAL MODELS FOR DE HAVILLAND CANADA "BUFFALO & TWIN OTTER" STOL TRANSPORTS

# VT: Vertical Tail
# HT: Horizontal Tail
# distance_HT : distance, wing quarter chord to horizontal tail quarter chord, ft
# distance_VT : distance, wing quarter chord to vertical tail quarter chord, ft
# Z_w : vertical distance of root quarter chord of wing below fuselage center line, ft
# Z_VT : vertival distance, aircraft CG to mean chord of vertival tail, ft
# Cl_beta_w_by_r1 : page 19 in the ref
# Cl_beta_w_by_r2 : page 19 in the ref
# y1 : distance, aircraft centerline to inboard end of ailereon, ft
# y2 : distance, aircraft centerline to outboard end of aileron, ft
# theta_0: angle of rotation, at equilibrium flight condition, of X_A axis from horizontal plane, rad or deg

# Secondary parameters
downwash_rate_alpha : (denote: dϵ/dα) rate of change of downwash angle at tail with change in angle of attack, rad/rad 
'''
import numpy as np


def add_VT_surface(label,CL_alpha_VT, ):
    
    pass


def add_HT_surface(label,CL_alpha_HT,distance_HT,area_HT,effi_factor_HT=1):

    pass

def add_wing():
    pass
class dynamic_stability_derivatives():
    '''
    Calculate maximul lift coefficient and plot lift curve using Raymer method
    '''
    def __init__(self,):

        self.CD_w
        self.CL_0
        self.CL_alpha_wing
        self.mac_wing
        self.area_main_wing
        self.wing_span
        self.tape_ratio_wing
        self.aspect_ratio_wing
        self.sweep_angle

        self.CL_alpha_HT
        self.distance_HT
        self.area_HT
        self.effi_factor_HT = 1
        
        self.CL_alpha_VT
        self.distance_VT
        self.area_VT
        self.tape_ratio_VT
        self.span_VT
        self.Z_VT
        self.effi_factor_VT = 1        
      
        # Secondary parameters
        self.horizontal_tal_vol = self.calc_horizontal_tail_vol(self.effi_factor_HT, self.distance_HT, self.area_HT, self.mac_wing, self.area_main_wing)
        self.downwash_rate_alpha = self.calc_downwash_rate_alpha( self.CL_alpha_wing, self.aspect_ratio_wing)

        # Calculation of aerodynamic coefficients 
        self.cx_alpha_dot = 0
        self.cx_q = 0

        self.cz_alpha_dot = self.calc_cz_alpha_dot(self.CL_alpha_HT, self.downwash_rate_alpha, self.horizontal_tal_vol)
        self.cz_q = self.calc_cz_q(self.CL_alpha_HT,self.horizontal_tal_vol)

        self.cm_alpha_dot= self.calc_cm_alpha_dot(self.distance_HT, self.mac_wing, self.cz_alpha_dot)
        self.cm_q = self.calc_cm_q( self.distance_HT, self.mac_wing, self.cz_q)

        self.cy_r = self.calc_cy_r(self.CL_alpha_VT, self.distance_VT, self.area_VT, self.wing_span, self.area_main_wing)
        self.cy_p = self.calc_cy_p(self.CL_0, self.sweep_angle)

        self.cl_r = self.calc_cl_r(CL_0, Z_VT, wing_span,cy_r)
        self.cl_p = self.calc_cl_p(CL_alpha_wing, self.tape_ratio_wing)

        self.cn_p = self.calc_cn_p(self.CL_0, self.CL_alpha_wing, self.aspect_ratio_wing, self.wing_span, self.distance_VT, self.cy_p)
        self.cn_r = self.calc_cn_r(self.distance_VT, self.wing_span, self.CD_w, self.cy_r)
        
    ###### Dynamic Stability Coefficients Estimating ########################
    # 1 CX due to alphadot
    def calc_cx_alpha_dot():
        return 0

    # 2 CX due to alphadot due to pithching rate
    def calc_cx_q():
        return 0

    # 3 CZ due to  the rate of  AOA ( tail contribution)
    def calc_cz_alpha_dot(CL_alpha_HT, downwash_rate_alpha,horizontal_tal_vol):
        return -2 * CL_alpha_HT * horizontal_tal_vol*downwash_rate_alpha

    # 4 CZ due  pithching rate
    def calc_cz_q(CL_alpha_HT,horizontal_tal_vol): #
        return -2 *CL_alpha_HT * horizontal_tal_vol

    # 5 CM due to  the rate of  AOA
    def calc_cm_alpha_dot(distance_HT, mac_wing,cz_alpha_dot): 
        return cz_alpha_dot *distance_HT/ mac_wing

    # 6 CM  due to pithching rate
    def calc_cm_q( distance_HT, mac_wing,cz_q):
        return cz_q*distance_HT/ mac_wing 

    # 7 CY due to yawing  rate
    def calc_cy_r(CL_alpha_VT, distance_VT, area_VT, wing_span, area_main_wing):
            return 2 * CL_alpha_VT * (distance_VT/wing_span) * (area_VT/area_main_wing)

    # 8 CY due to rolling rate
    # def calc_cy_p(CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT):
    #     return -2/3 * CL_alpha_VT * (span_VT/wing_span) * (area_VT/area_main_wing) * ((1+ 2*tape_ratio_VT)/(1+tape_ratio_VT)) 
    
    def calc_cy_p(CL_0,sweep_angle):
        return CL_0*(aspect_ratio_wing + np.cos(sweep_angle))/(aspect_ratio_wing + 4*np.cos(sweep_angle))*np.tan(sweep_angle)

    # 9 C_roll due to yawing rate  
    def calc_cl_r (CL_0, Z_VT, wing_span,cy_r):
        return CL_0 / 4 + cy_r* Z_VT/wing_span 

    # 10 C_roll  due to rolling rate
    def calc_cl_p (CL_alpha_wing,tape_ratio_wing):
        return -CL_alpha_wing/12 * (1 + 3*tape_ratio_wing)/(1 + tape_ratio_wing) 

    # 11 Cn due to rolling rate
    
    def calc_cn_p(CL_0,CL_alpha_wing, aspect_ratio_wing,wing_span,distance_VT,cy_p):
        return -CL_0/4 * (1 - CL_alpha_wing/(np.pi*aspect_ratio_wing)) - cy_p * distance_VT /wing_span

    # 12 Cn due to yawing rate
    def calc_cn_r(distance_VT, wing_span, CD_w,cy_r):
        return -(distance_VT/wing_span)*cy_r - CD_w/4 

    #  Other secondary parameters 

    #  HT volume ratio
    def calc_horizontal_tail_vol(effi_factor_HT,distance_HT, area_HT, mac_wing, area_main_wing):
        return effi_factor_HT * (distance_HT / mac_wing) * (area_HT / area_main_wing)

    def calc_downwash_rate_alpha(CL_alpha_wing,aspect_ratio_wing): # 
        return 2*CL_alpha_wing/(3.14*aspect_ratio_wing)

    ######################################### Static Stability Coefficients Estimating #######################################################

    # def calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT):
    #     return CL_alpha_HT * elevator_effectiveness * (area_HT / area_main_wing) *effi_factor_HT

    # def calc_cm_delta_e(distance_HT, mac_wing,CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT):
    #     return calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT) *distance_HT/ mac_wing

    # def calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness):
    #     return -CL_alpha_VT * (area_VT/area_main_wing) * rudder_effectiveness

    # def calc_cn_sigma_a():
    #     return 0

    # def calc_cn_sigma_r(distance_VT,wing_span,CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness):
    #     return -self.calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness) * distance_VT/wing_span

    # def calc_cl_beta (self, Cl_beta_w_by_r1, wing_dihedral_angle, Cl_beta_w_by_r2,area_VT, area_main_wing, 
    #             CL_alpha_VT, downwash_rate_beta,aspect_ratio_wing,Z_w,wing_span,fuselage_height,fuselage_width,Z_VT):
    #     return ((Cl_beta_w_by_r1 - Cl_beta_w_by_r2) * wing_dihedral_angle + 1.2*np.sqrt(aspect_ratio_wing) * (Z_w/wing_span) 
    #             * (fuselage_height + fuselage_width)/wing_span + self.calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta)*Z_VT/wing_span)

    # def calc_cl_sigma_a(CL_alpha_wing, aileron_effectiveness, area_main_wing, wing_span, y1, y2):
    #     return (2 * CL_alpha_wing * aileron_effectiveness / (area_main_wing * wing_span)) * (7.41*(y2**2 - y1**2) - 0.062 * (47**3 - 33**3))

    # def calc_cl_sigma_r(self, CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness, wing_span, Z_VT):
    #     return self.calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness) * Z_VT/wing_span

    # def calc_cx_u(CD_0,CL_0,theta_0 = 0):
    #     return - 3 * CD_0 - CL_0 * np.tan(theta_0) 

    # def calc_cx_alpha(CL_0, CL_alpha_wing, aspect_ratio_wing, aircraft_effi_factor):
    #     return CL_0 - 2* CL_0 * CL_alpha_wing / (np.pi * aircraft_effi_factor * aspect_ratio_wing)

    # def calc_cz_u(CL_0):
    #     return -2 * CL_0

    # def calc_cz_alpha(CL_alpha_wing, CD_0):
    #     return -(CL_alpha_wing + CD_0)

    # def calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta):
    #     return -CL_alpha_VT * area_VT / area_main_wing * (1 - downwash_rate_beta)

    # def calc_cn_beta(self,distance_VT, wing_span, vol_fuse,area_main_wing,fuselage_height,fuselage_width,area_VT, CL_alpha_VT, downwash_rate_beta):
    #     return (-self.calc_cy_beta(area_VT, area_main_wing, CL_alpha_VT, downwash_rate_beta) * distance_VT/wing_span 
    #             - 1.3 * vol_fuse/(area_main_wing*wing_span) * fuselage_height/fuselage_width)




if __name__ == "__main__":

    effi_factor_HT = 1.0
    distance_HT= 46
    area_HT = 233
    mac_wing = 10.1
    area_main_wing = 945
    CL_alpha_HT = 3.5
    CL_alpha_VT = 2.5
    span_VT = 14
    wing_span = 96
    area_VT =152
    tape_ratio_VT = 0.64
    distance_VT = 44
    CD_w = 0.025
    CL_0 = 0.77
    CD_0 = 0.057
    aspect_ratio_wing = 9.75
    Z_VT = 10
    Z_w = -3.5
    y1 = 33
    y2 = 47
    cl_p = -0.53
    CL_alpha_wing = 5.2
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



    cz_q = calc_cz_q(CL_alpha_HT, effi_factor_HT, distance_HT, area_HT, mac_wing, area_main_wing)
    cm_q = calc_cm_q(distance_HT, mac_wing,CL_alpha_HT, effi_factor_HT, area_HT, area_main_wing)
    cy_p = calc_cy_p(CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT)
    cy_r = calc_cy_r(CL_alpha_VT, distance_VT, area_VT, wing_span, area_main_wing)
    cn_p = calc_cn_p(CL_0,CL_alpha_wing, aspect_ratio_wing,CL_alpha_VT,span_VT,wing_span,area_VT,area_main_wing,tape_ratio_VT,distance_VT)
    cn_r = calc_cn_r(distance_VT, wing_span, CD_w, CL_alpha_VT, area_VT, area_main_wing)
    cl_p_value = calc_cl_p(cl_p)
    cl_r = calc_cl_r(CL_0, Z_VT, wing_span,CL_alpha_VT,distance_VT,area_VT,area_main_wing)
    cz_alpha_dot = calc_cz_alpha_dot(CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, distance_HT, area_HT, mac_wing, area_main_wing)
    cz_delta_e = calc_cz_delta_e(CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT)
    cm_alpha_dot = calc_cm_alpha_dot(distance_HT, mac_wing, CL_alpha_HT, downwash_rate_alpha, effi_factor_HT, area_HT, area_main_wing)
    cm_delta_e = calc_cm_delta_e(distance_HT, mac_wing,CL_alpha_HT, elevator_effectiveness, effi_factor_HT, area_main_wing, area_HT)
    cy_sigma_r = calc_cy_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness)
    cn_sigma_a = calc_cn_sigma_a()
    cn_sigma_r = calc_cn_sigma_r(distance_VT,wing_span,CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness)
    cl_sigma_a = calc_cl_sigma_a(CL_alpha_wing, aileron_effectiveness, area_main_wing, wing_span, y1, y2)
    cl_sigma_r = calc_cl_sigma_r(CL_alpha_VT, area_main_wing, area_VT, rudder_effectiveness, wing_span, Z_VT)
    cx_u = calc_cx_u(CD_0, CL_0, theta_0)
    cx_alpha = calc_cx_alpha(CL_0,CL_alpha_wing, aspect_ratio_wing, aircraft_effi_factor)
    cz_u = calc_cz_u(CL_0)
    cz_alpha = calc_cz_alpha(CL_alpha_wing,CD_0)
    cy_beta = calc_cy_beta(area_VT,area_main_wing,CL_alpha_VT,downwash_rate_beta)
    cn_beta = calc_cn_beta(distance_VT,wing_span,vol_fuse,area_main_wing,fuselage_height,fuselage_width,area_VT,CL_alpha_VT,downwash_rate_beta) 
    cl_beta = calc_cl_beta (Cl_beta_w_by_r1, wing_dihedral_angle, Cl_beta_w_by_r2,area_VT, area_main_wing, 
                CL_alpha_VT, downwash_rate_beta,aspect_ratio_wing,Z_w,wing_span,fuselage_height,fuselage_width,Z_VT)

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










        
        
    
    
    
    
        