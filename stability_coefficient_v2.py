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
CD_w : wing drag coef.

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
Z_VT : vertical distance, aircraft CG to mean chord of vertival tail, ft
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

#---------------------------------WING----------------------------
class main_wing():
    '''
    Calculate contribution of single main_wing
    '''
    def __init__(self,CL_0, CL_alpha_wing, CD_w, sweep_angle, aspect_ratio_wing,tape_ratio_wing):
        self.cy_p = self.calc_cy_p_wing(CL_0,sweep_angle, aspect_ratio_wing)
        self.cl_r = self.calc_cl_r_wing (CL_0)
        self.cl_p = self.calc_cl_p_wing (CL_alpha_wing,tape_ratio_wing)
        self.cn_p= self.calc_cn_p_wing(CL_0,CL_alpha_wing, aspect_ratio_wing)
        self.cn_r = self.calc_cn_r_wing(CD_w)

    #CY due to rolling rate - - wing contribution
    def calc_cy_p_wing(self,CL_0,sweep_angle, aspect_ratio_wing):
        return CL_0*(aspect_ratio_wing + np.cos(sweep_angle))/(aspect_ratio_wing + 4*np.cos(sweep_angle))*np.tan(sweep_angle)
        
    # 10 C_roll  due to rolling rate- wing contribution
    def calc_cl_p_wing (self,CL_alpha_wing,tape_ratio_wing):
        return -CL_alpha_wing/12 * (1 + 3*tape_ratio_wing)/(1 + tape_ratio_wing) 

    def calc_cl_r_wing (self,CL_0):
        return CL_0 / 4  

    # 11 Cn due to rolling rate
    def calc_cn_p_wing(self,CL_0,CL_alpha_wing, aspect_ratio_wing):
        return -CL_0/4 * (1 - CL_alpha_wing/(np.pi*aspect_ratio_wing))     

    # 12 Cn due to yawing rate
    def calc_cn_r_wing(self,CD_w):
        return  - CD_w/4     

class HT_surface():
    '''
    Calculate contribution of single HT_surface
    '''
    def __init__(self,CL_alpha_HT, distance_HT, area_HT, CL_alpha_wing, mac_wing, area_wing, aspect_ratio_wing, effi_factor_HT = 1):
#---------------------------------HT tail----------------------------
        self.horizontal_tal_vol = self.calc_horizontal_tail_vol(effi_factor_HT, distance_HT, area_HT, mac_wing, area_wing)
        self.downwash_rate_alpha = self.calc_downwash_rate_alpha(CL_alpha_wing, aspect_ratio_wing)   
    
        self.cz_alpha_dot = self.calc_cz_alpha_dot(CL_alpha_HT, self.downwash_rate_alpha, self.horizontal_tal_vol)
        self.cz_q = self.calc_cz_q(CL_alpha_HT, self.horizontal_tal_vol)
        self.cm_alpha_dot= self.calc_cm_alpha_dot(distance_HT, mac_wing, self.cz_alpha_dot)
        self.cm_q = self.calc_cm_q( distance_HT, mac_wing, self.cz_q)

    #  HT volume ratio
    def calc_horizontal_tail_vol(self, effi_factor_HT, distance_HT, area_HT, mac_wing, area_wing):
        return effi_factor_HT * (distance_HT / mac_wing) * (area_HT / area_wing)  

    def calc_downwash_rate_alpha(self, CL_alpha_wing,aspect_ratio_wing): # 
        return 2*CL_alpha_wing/(3.14*aspect_ratio_wing)
        
    # 3 CZ due to  the rate of  AOA ( tail contribution)
    def calc_cz_alpha_dot(self, CL_alpha_HT, downwash_rate_alpha, horizontal_tal_vol):
        return -2 * CL_alpha_HT * horizontal_tal_vol*downwash_rate_alpha

    # 4 CZ due  pithching rate - tail contribution
    def calc_cz_q(self, CL_alpha_HT, horizontal_tal_vol): #
        return -2 *CL_alpha_HT * horizontal_tal_vol

    # 5 CM due to  the rate of  AOA - tail contribution
    def calc_cm_alpha_dot(self, distance_HT, mac_wing, cz_alpha_dot): 
        return cz_alpha_dot *distance_HT/ mac_wing

    # 6 CM  due to pithching rate - tail contribution
    def calc_cm_q(self, distance_HT, mac_wing,cz_q):
        return cz_q*distance_HT/ mac_wing     


class VT_surface():
    '''
    Calculate contribution of single VT_surface
    Z_wing :  distance, parallel to z axis, from wing root quarter chord to fuselage centerline
    '''
    def __init__(self,CL_alpha_VT, distance_VT, area_VT, span_VT, tape_ratio_VT, Z_VT, CL_0, wing_span, area_wing,aspect_ratio_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth):
#---------------------------------VT tail----------------------------
        self.sidewash_rate_beta = self.calc_sidewash_rate_beta(area_VT,area_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth, aspect_ratio_wing )
        self.cy_beta = self.calc_cy_beta_VT_tail(CL_alpha_VT, area_VT, area_wing , self.sidewash_rate_beta)

        self.cy_r = self.calc_cy_r_VT_tail(distance_VT, wing_span, self.cy_beta)
        self.cy_p = self.calc_cy_p_VT_tail(CL_alpha_VT, span_VT, wing_span,area_VT, area_wing, tape_ratio_VT)
        self.cl_r = self.calc_cl_r(CL_0, Z_VT, wing_span, self.cy_r)
        self.cn_p = self.calc_cn_p_VT_tail(wing_span, distance_VT, self.cy_p)
        self.cn_r = self.calc_cn_r_VT_tail(distance_VT, wing_span, self.cy_r)

    def calc_sidewash_rate_beta (self, area_VT, area_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth, aspect_ratio_wing):
        return 0.724 + 300.6*(area_VT/area_wing)/(1 + np.cos(sweep_at_quater_mac_wing))+ 0.4*Z_wing/fulselage_depth + 0.009*aspect_ratio_wing
    
    def calc_vertical_tail_vol(self, distance_VT, area_VT, mac_wing, area_wing):
        return  (distance_VT/mac_wing) * (area_VT / area_wing)      

    def calc_cy_beta_VT_tail (self, CL_alpha_VT, area_VT, area_wing , sidewash_rate_beta):
        return - CL_alpha_VT*area_VT/area_wing *sidewash_rate_beta

    # 7 CY due to yawing  rate -tail contribution
    def calc_cy_r_VT_tail (self, distance_VT, wing_span, cy_beta):
            return -2 *(distance_VT/wing_span) * cy_beta

    # 8 CY due to rolling rate - - tail contribution
    def calc_cy_p_VT_tail (self, CL_alpha_VT, span_VT, wing_span, area_VT, area_wing, tape_ratio_VT):
        return -2/3 * CL_alpha_VT * (span_VT/wing_span) * (area_VT/area_wing) * ((1+ 2*tape_ratio_VT)/(1+tape_ratio_VT))        

    # 9 C_roll due to yawing rate  
    def calc_cl_r (self, CL_0, Z_VT, wing_span, cy_r):
        return  CL_0/4 + cy_r* Z_VT/wing_span 

    # 11 Cn due to rolling rate
    def calc_cn_p_VT_tail(self, wing_span, distance_VT, cy_p):
        return - cy_p * distance_VT /wing_span      

    # Cn due to yowing rate
    def calc_cn_r_VT_tail(self, distance_VT, wing_span, cy_r):
        return -(distance_VT/wing_span)*cy_r 
    


if __name__ == "__main__":
    #Kp2 configurtion
    CL_0                        = 0.2960
    CL_alpha_wing               = 4.3544
    CD_w                        = 0.0045
    sweep_angle                 = 20*np.pi/180 #radian, at leading edge
    aspect_ratio_wing           = 6.4418
    mac_wing                    = 0.2988 #m
    wing_span                   = 1.9248 #m
    area_wing                   = 0.534041 #m^2
    sweep_at_quater_mac_wing    = 16.54*np.pi/180 #radian
    Z_wing                      = -0.1371
    fulselage_depth             = 0.3304
    tape_ratio_wing             = 2.8512
    # main wing

    CL_alpha_HT_1               = 2.5456
    distance_HT_1               = 0.9969 #m
    area_HT_1                   = 0.134241 #m^2

    # slope angle of HT 2 & 3 = 50 degree
    CL_alpha_HT_2               = 5.7547
    distance_HT_2               = 0.9130
    area_HT_2                   = 0.059720 #m^2

    CL_alpha_HT_3               = 5.7547
    distance_HT_3               = 0.9130
    area_HT_3                   = 0.059720 # m^2
    
    CL_alpha_HT_4               = 2.5456
    distance_HT_4               = 0.9969 #m
    area_HT_4                   = 0.134241 #m^2

    CL_alpha_VT_1               = 4.7956
    distance_VT_1               = 0.5668 #m
    area_VT_1                   = 0.03941 #m^2
    span_VT_1                   = 0.2346 #m
    tape_ratio_VT_1             = 0.7998
    Z_VT_1                      = 0.3591

    CL_alpha_VT_2               = 5.7547
    distance_VT_2               = 0.9131
    area_VT_2                   = 0.07092 #m^2
    span_VT_2                   = 0.2869 #m
    tape_ratio_VT_2             = 0.9180
    Z_VT_2                      = 0.3860

    CL_alpha_VT_3               = 5.7547
    distance_VT_3               = 0.9131
    area_VT_3                   = 0.07092 #m^2
    span_VT_3                   = 0.2869 #m
    tape_ratio_VT_3             = 0.9180
    Z_VT_3                      = 0.3860

    CL_alpha_VT_4               = 4.7956
    distance_VT_4               = 0.5668 #m
    area_VT_4                   = 0.03941 #m^2
    span_VT_4                   = 0.2346
    tape_ratio_VT_4             = 0.7998
    Z_VT_4                      = 0.3591

    wing = main_wing(CL_0, CL_alpha_wing, CD_w, sweep_angle, aspect_ratio_wing, tape_ratio_wing)

    # rearwing- HT tail 1
    h_tail_1 = HT_surface(CL_alpha_HT_1, distance_HT_1, area_HT_1, CL_alpha_wing, mac_wing, area_wing, aspect_ratio_wing, effi_factor_HT = 1)
    # V-tail 
        # HT 2
    h_tail_2 = HT_surface(CL_alpha_HT_2, distance_HT_2, area_HT_2, CL_alpha_wing, mac_wing, area_wing, aspect_ratio_wing, effi_factor_HT = 1)        
        # HT 3
    h_tail_3 = HT_surface(CL_alpha_HT_3, distance_HT_3, area_HT_3, CL_alpha_wing, mac_wing, area_wing, aspect_ratio_wing, effi_factor_HT = 1)  
        # VT 1
    v_tail_1 = VT_surface(CL_alpha_VT_1, distance_VT_1, area_VT_1, span_VT_1, tape_ratio_VT_1, Z_VT_1, CL_0, wing_span, area_wing,aspect_ratio_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth)
        # VT 2
    v_tail_2 = VT_surface(CL_alpha_VT_2, distance_VT_2, area_VT_2, span_VT_2, tape_ratio_VT_2, Z_VT_2, CL_0, wing_span, area_wing,aspect_ratio_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth)
        # VT 2  
    # Winglet
        # VT 3
    v_tail_3 = VT_surface(CL_alpha_VT_3, distance_VT_3, area_VT_3, span_VT_3, tape_ratio_VT_3, Z_VT_3, CL_0, wing_span, area_wing, aspect_ratio_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth)
        # VT 4
    v_tail_4 = VT_surface(CL_alpha_VT_4, distance_VT_4, area_VT_4, span_VT_4, tape_ratio_VT_4, Z_VT_4, CL_0, wing_span, area_wing, aspect_ratio_wing, sweep_at_quater_mac_wing, Z_wing, fulselage_depth)

    ## result

            # Calculation of aerodynamic coefficients 
    cx_alpha_dot = 0
    cx_q = 0

    cz_alpha_dot = h_tail_1.cz_alpha_dot + h_tail_2.cz_alpha_dot + h_tail_3.cz_alpha_dot
    cz_q = h_tail_1.cz_q + h_tail_2.cz_q + h_tail_3.cz_q

    cm_alpha_dot=  h_tail_1.cm_alpha_dot + h_tail_2.cm_alpha_dot + h_tail_3.cm_alpha_dot
    cm_q = h_tail_1.cm_q + h_tail_2.cm_q + h_tail_3.cm_q

    cy_r = v_tail_1.cy_r + v_tail_2.cy_r + v_tail_3.cy_r + v_tail_4.cy_r 
    cy_p = v_tail_1.cy_p + v_tail_2.cy_p + v_tail_3.cy_p + v_tail_4.cy_p + wing.cy_p

    cl_r = v_tail_1.cy_r + v_tail_2.cy_r + v_tail_3.cy_r + v_tail_4.cy_r + wing.cl_r
    cl_p =  wing.cl_p

    cn_p = v_tail_1.cn_p + v_tail_2.cn_p + v_tail_3.cn_p + v_tail_4.cn_p + wing.cn_p
    cn_r = v_tail_1.cn_r + v_tail_2.cn_r + v_tail_3.cn_r + v_tail_4.cn_r + wing.cn_r
    
    print('\ncx_alpha_dot =',cx_alpha_dot)
    print('cx_q =',cx_q)
    print('cz_alpha_dot =',cz_alpha_dot)
    print('cz_q =',cz_q)
    print('cm_q =',cm_q)
    print('cy_r =',cy_r)
    print('cy_p =',cy_p)
    print('cl_r =',cl_r)
    print('cl_p =',cl_p)
    print('cn_p =',cn_p)
    print('cn_r =',cn_r)
