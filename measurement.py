# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 22:33:18 2023

@author: nghianguyen
"""
import numpy as np

class CONFIG_KP2C:
    def __init__(self):
        self.main_wing_sweep_angle_quater_chord = None #degree
        self.main_wing_aspect_ratio             = None 
        self.main_wing_mac                      = None #m
        self.main_wing_span                     = None #m
        self.main_wing_area                     = None #m
        self.main_wing_z                        = None #m
        self.fuselage_depth                     = None #m
        self.ht_1_area                          = None #m^2
        self.ht_2_area                          = None #m^2
        self.ht_3_area                          = None #m^2
        self.ht_4_area                          = None #m^2
        self.vt_1_area                          = None #m^2
        self.vt_2_area                          = None #m^2
        self.vt_3_area                          = None #m^2
        self.vt_4_area                          = None #m^2
        self.vt_1_span                          = None #m
        self.vt_2_span                          = None #m
        self.vt_3_span                          = None #m
        self.vt_4_span                          = None #m
        self.vt_1_taper_ratio                    = None
        self.vt_2_taper_ratio                    = None
        self.vt_3_taper_ratio                    = None
        self.vt_4_taper_ratio                    = None
        
        
    def calc_distance_ht_1(self,x_mainwing_quarter_chord,x_ht_1_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_ht_1_quarter_chord))
    
    def calc_distance_ht_2(self,x_mainwing_quarter_chord,x_ht_2_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_ht_2_quarter_chord))
    
    def calc_distance_ht_3(self,x_mainwing_quarter_chord,x_ht_3_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_ht_3_quarter_chord))
    
    def calc_distance_ht_4(self,x_mainwing_quarter_chord,x_ht_4_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_ht_4_quarter_chord))
    
    def calc_distance_vt_1(self,x_mainwing_quarter_chord,x_vt_1_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_vt_1_quarter_chord))
    
    def calc_distance_vt_2(self,x_mainwing_quarter_chord,x_vt_2_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_vt_2_quarter_chord))
    
    def calc_distance_vt_3(self,x_mainwing_quarter_chord,x_vt_3_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_vt_3_quarter_chord))
    
    def calc_distance_vt_4(self,x_mainwing_quarter_chord,x_vt_4_quarter_chord):
        return np.sqrt(np.abs(x_mainwing_quarter_chord -x_vt_4_quarter_chord))
    
    def calc_z_vt_1(self,z_CG, z_vt_1_mean_chord):
        return np.sqrt(np.abs(z_CG[2] -z_vt_1_mean_chord))
    
    def calc_z_vt_2(self,z_CG, z_vt_2_mean_chord):
        return np.sqrt(np.abs(z_CG[2] -z_vt_2_mean_chord))
    
    def calc_z_vt_3(self,z_CG, z_vt_3_mean_chord):
        return np.sqrt(np.abs(z_CG[2] -z_vt_3_mean_chord))
    
    def calc_z_vt_4(self,z_CG, z_vt_mean_4_chord):
        return np.sqrt(np.abs(z_CG[2] -z_vt_mean_4_chord))
    
if __name__ == "__main__":
    z_CG                        = [1.3605, 0, 1.0283] 
    x_mainwing_quarter_chord    = 1.1317
    x_ht_1_quarter_chord        = 2.1256
    x_ht_2_quarter_chord        = 1.9654
    x_ht_3_quarter_chord        = 1.9654
    x_ht_4_quarter_chord        = 2.1256
    x_vt_1_quarter_chord        = 1.4530
    x_vt_2_quarter_chord        = 1.9655
    x_vt_3_quarter_chord        = 1.9655
    x_vt_4_quarter_chord        = 1.4530
    z_vt_1_mean_chord           = 1.1573
    z_vt_2_mean_chord           = 1.1773
    z_vt_3_mean_chord           = 1.1773
    z_vt_4_mean_chord           = 1.1573
    
    para = CONFIG_KP2C()
    para.main_wing_sweep_angle_quater_chord = 16.405 
    para.main_wing_area                     = 53.4041
    para.main_wing_span                     = 1.9248
    para.main_wing_mac                      = 0.2988
    para.main_wing_aspect_ratio             = 6.4418
    para.main_wing_z                        = -0.1371
    para.fuselage_depth                     = 0.3304
    para.ht_1_area                          = 13.4241
    para.ht_2_area                          = 5.9720
    para.ht_3_area                          = 5.9720
    para.ht_4_area                          = 13.4241
    para.vt_1_area                          = 3.941
    para.vt_2_area                          = 7.092
    para.vt_3_area                          = 7.092
    para.vt_4_area                          = 3.941
    para.vt_1_span                          = 0.2346
    para.vt_2_span                          = 0.2869
    para.vt_3_span                          = 0.2869
    para.vt_4_span                          = 0.2346
    para.vt_1_taper_ratio                   = 0.7998
    para.vt_2_taper_ratio                   = 0.9180
    para.vt_3_taper_ratio                   = 0.9180
    para.vt_4_taper_ratio                   = 0.7998
    
    #### calc required parameters ###
    
    # distance, wing quarter chord to horizontal tail quarter chord, m
    lt_1 = para.calc_distance_ht_1(x_mainwing_quarter_chord,x_ht_1_quarter_chord)
    lt_2 = para.calc_distance_ht_2(x_mainwing_quarter_chord,x_ht_2_quarter_chord)
    lt_3 = para.calc_distance_ht_3(x_mainwing_quarter_chord,x_ht_3_quarter_chord)
    lt_4 = para.calc_distance_ht_4(x_mainwing_quarter_chord,x_ht_4_quarter_chord)
    
    # distance, wing quarter chord to vertical tail quarter chord, m
    lf_1 = para.calc_distance_vt_1(x_mainwing_quarter_chord, x_vt_1_quarter_chord)
    lf_2 = para.calc_distance_vt_2(x_mainwing_quarter_chord, x_vt_2_quarter_chord)
    lf_3 = para.calc_distance_vt_3(x_mainwing_quarter_chord, x_vt_3_quarter_chord)
    lf_4 = para.calc_distance_vt_4(x_mainwing_quarter_chord, x_vt_4_quarter_chord)
    
    # vertical distance, aircraft CG to mean chord of the vertical tail, m
    z1 = para.calc_z_vt_1(z_CG, z_vt_1_mean_chord)
    z2 = para.calc_z_vt_2(z_CG, z_vt_2_mean_chord)
    z3 = para.calc_z_vt_3(z_CG, z_vt_3_mean_chord)
    z4 = para.calc_z_vt_4(z_CG, z_vt_4_mean_chord)
    
    
    
    
    
    