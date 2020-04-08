import demo_data
import demo_data_classical_Thellier
import demo_data_MW
import demo_data_MSP
import demo_data_pTh_GF
import demo_data_pTh_NAA
import graphing
import graphing_pTh
import helpers
import helpers_curv
import sc_main
import sc_config

import math

import MSP_param
import MSP_bootstrap


def pp_dict(d, indent=0):
    for key, value in d.items():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pp_dict(value, indent + 1)
        else:
            print('\t' * (indent + 1) + str(value))


#############
# DEMO DATA #
#############
# print(demo_data.get("RE04"))
# print(demo_data.get("RE16"))
# print(demo_data.get("RE16", True))
# print(demo_data.get("ET187", True))
# print(demo_data.get("garbage"))
# print(demo_data.get("garbage", True))

# print(demo_data_pTh_GF.get("RE01", True))


############
# PLOTTING #
###########
# specimen = demo_data.get("RE04")

specimen = demo_data_classical_Thellier.get("QCR135A")

print(graphing.plot_arai(specimen))
# print(graphing.plot_zijderveld(specimen))
# print(graphing.plot_equalarea(specimen))
# print(graphing.plot_magnitisation(specimen))


#############
#  HELPERS  #
#############
# print(helpers.cart2dir(2,2,2))
# print(helpers.dir2cart(2,2,2))
# print(helpers.get_angle_diff([1,1,1], [-1,-1,-1]))

##############
# pTh params #
##############
# specimen = demo_data_pTh_NAA.get("RE01")

# print(pTh_param.B12ARM_param_pth_NAA(specimen))

# print(graphing_pTh.plot_arai_pTh(specimen))
# print(graphing_pTh.plot_zijderveld_pTh(specimen))
# print(graphing_pTh.plot_equalarea_pTh(specimen))
# print(graphing_pTh.plot_magnitisation_pTh(specimen))
# print(graphing_pTh.B12ARM_param_pth(specimen))
# print(graphing_pTh.plot_pTh_graphs(specimen))

##############
# MSP params #
##############
# site = demo_data_MSP.get("PI320_db")
# alpha_MSP = 0.5

# print(MSP_param.MSP_params_calc(site, alpha_MSP))


#################
# MSP Bootstrap #
################
# site_MSP = demo_data_MSP.get("PI320_db_only")       # change to "PI320_db" to test with the None type thingy
# selection_MSP = demo_data_MSP.get("PI320_db_only", True)
# alpha_MSP = 0.5
# NumCycles_MSP = 1000
# cutOff_MSP = 95

# print(MSP_bootstrap.MSP_boots_calc(site_MSP, selection_MSP, alpha_MSP, NumCycles_MSP, cutOff_MSP))


######################
# SELECTION CRITERIA #
# ####################
specimen_thermal = demo_data.get("BJY2")
selection_thermal = demo_data.get("BJY2", True)

specimen_TT = demo_data_classical_Thellier.get("QCR135A")
selection_TT = demo_data_classical_Thellier.get("QCR135A", True)

specimen_MW = demo_data_MW.get("MNSAMS106B")
selection_MW = demo_data_MW.get("MNSAMS106B", True)

specimen_pTh_GF = demo_data_pTh_GF.get("RE01")
selection_pTh_GF = demo_data_pTh_GF.get("RE01", True)

specimen_pTh_NAA = demo_data_pTh_NAA.get("RE01")
selection_pTh_NAA = demo_data_pTh_NAA.get("RE01", True)

# suitcase1 = sc_main.thermal_thellier(specimen_thermal, selection_thermal)
# suitcase2 = sc_main.microwave_thellier(specimen_MW, selection_MW)
# suitcase5 = sc_main.pseudo_thellier(specimen_pTh_NAA, selection_pTh_NAA)  # works now also for GF and NAA format (specimen_pTh_GF, selection_pTh_GF) (specimen_pTh_NAA, selection_pTh_NAA)

# suitcase to test Thellier-Thellier, only for testing purposes, uses the same suitcase as Thermal Thellier
suitcase4 = sc_main.thermal_thellier(specimen_TT, selection_TT)


# # MSP #####

# site_MSP = demo_data_MSP.get("PI320_db_only")       # change to "PI320_db" to test with the None type thingy
# selection_MSP = demo_data_MSP.get("PI320_db_only", True)
# alpha_MSP = 0.5
# NumCycles_MSP = 1000
# cutOff_MSP = 95

# suitcase3 = sc_main.MSP(site_MSP, selection_MSP, alpha_MSP, NumCycles_MSP, cutOff_MSP)


# ##### print #####
# pp_dict(suitcase1["preprocessed"]["basics"])#

pp_dict(suitcase4["anisotropy_statistics"])
