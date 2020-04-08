
# WARNING: this deepcopy only works with dictionaries!!


def deepcopy(elem):
    if (elem == None):
        return None
    else:
        new = {}
        for key, value in elem.items():
            new[key] = deepcopy(value)
        return new


def new_Thellier_suitcase():
    return deepcopy(Thellier_suitcase)


def new_TT_suitcase():
    return deepcopy(TT_suitcase)


def new_pseudo_Thellier_suitcase():
    return deepcopy(pseudo_Thellier_suitcase)


def new_MSP_suitcase():
    return deepcopy(MSP_suitcase)


# suitcase is the data structure that will be passed to all the
# steps in the selection criteria pipeline. It consists of
# multiple levels of dictionaries.
Thellier_suitcase = {
    "input": {
        "specimen": None,
        "specimen_selection": None
    },

    "preprocessed": {
        "msrmnts": {
            "zerofield_steps": None,
            "zerofield_steps_all": None,
            "infield_steps": None,
            "infield_steps_all": None,
            "infield_antiparallel_steps": None,
            "infield_antiparallel_steps_all": None,
            "ptrm_check_steps": None,
            "tail_check_steps": None,
            "add_check_steps": None,
        },

        "aniso_trm": {
            "orderd": None,
            "x+": None,
            "x-": None,
            "y+": None,
            "y-": None,
            "z+": None,
            "z-": None,
            "check": None,
        },

        "basics": {
            "nrm0": None,                   # = 1.04e-5 (one number only)
            "ptrm": None,                   # [{step, x, y, z}]
            "NRM_rem": None,                # [{step, x, y, z}]
            "NRM_vec_all": None,            # [x,y,z]
            "ptrm_vec_all": None,           # [x,y,z]
            "pstep_all": None,              # [step]
            "y_nrm": None,                  # [y]
            "xy_temp": None,                # [step]
            "NRM_vec_select": None,         # [x,y,z]
            "x_ptrm": None,                 # [x]
            "x_temp": None,                 # [step]
            "ptrm_gained_vec": None,        # [x,y,z]
            "xBar": None,                   # = 3.04e-6 (one number only)
            "yBar": None,                   # = 1.04e-6 (one number only)
            "x_ptrm_all": None,             # [x]
            "x_temp_all": None,             # [step]
            "y_nrm_all": None,              # [x]
            "y_temp_all": None,             # [step]
        },
        "field_basics": {
            "Blab": None,                   # = 50 (one number only)
            "field_dir_vec_initial": None,  # = x, y, z (one line only)
            "field_dir_vec": None,          # = x, y, z (one line only)
        },
        "checks": {
            "y_tail_check": None,           # [y]
            "x_tail_check": None,           # [x]
            "y_temp_tail_check": None,      # [step]
            "x_ptrm_check": None,           # [x]
            "y_ptrm_check": None,           # [y]
            "x_temp_ptrm_check": None,      # [step]
            "ptrm_check": None,             # [{step, x, y, z}]
            "add_check_step": None,         # [step]
            "add_check_vec": None,          # [x,y,z]
            "tail_check_vec": None,         # [x,y,z]
            "SCAT_ptrm_check_step": None,   # [step]
            "SCAT_tail_check_step": None,   # [step]
        }
    },

    "arai_statistics": {
        "PI_est": {
            "n": None,                      # = 8 (one number only)
            "b_slope": None,                # = 1.02 (one number only)
            "b_sigma": None,                # = 0.04 (one number only)
            "b_beta": None,                 # = 1.02 (one number only)
            "SE_anc": None,                 # = 1.02 (one number only)
            "sumx": None,                   # = 1.02 (one number only)
            "sumy": None,                   # = 1.02 (one number only)
            "sum_xy": None,                 # = 1.02 (one number only)
        },

        "PI_Banc_est": {
            "B_anc": None,                  # = 50.0 (one number only)
        },

        "intercept_stats": {
            "Yint": None,                   # = 6.02e-6 (one number only)
            "Xint": None                    # = 7.02e-6 (one number only)
        },
        "VDS_stats": {
            "VDS": None,                    # = 1.12e-5 (one number only)
            "FRAC": None,                   # = 0.42 (one number only)
            "GAP_MAX": None                 # = 0.20 (one number only)
        },
        "xy_prime_stat": {
            "x_prime": None,                # [x]
            "y_prime": None,                # [y]
        },
        "delta_prime_stat": {
            "delta_x_prime": None,          # = 5.258e-6 (one number only)
            "delta_y_prime": None,          # = 5.352e-6 (one number only)
        },
        "f_stat": {
            "f": None,                      # = 0.795 (one number only)
        },
        "f_VDS_stat": {
            "f_VDS": None,                  # = 0.475 (one number only)
        },
        "g_stat": {
            "g": None,                      # = 0.857 (one number only)
            "g_lim": None,                  # = 0.818 (one number only)
        },
        "q_stat": {
            "q": None,                      # = 13.7 (one number only)
        },
        "w_stat": {
            "w": None,                      # = (one number only)
        },
        "rsq_stat": {
            "rsq_corr": None,               # = (one number only)
            "rsq_det": None,                # = (one number only)
        },
    },

    "directional_statistics": {
        "mean_dir_stat": {
            "Mdec_free": None,              # = 88.81  (one number only)
            "Minc_free": None,              # = -51.33 (one number only)
            "MAD_free": None,               # = 6.26   (one number only)
            "Mdec_anc": None,               # = 96.01  (one number only)
            "Minc_anc": None,               # = -50.829(one number only)
            "MAD_anc": None,                # = 3.483  (one number only)
            "x": None,                      # [x]
            "y": None,                      # [y]
            "z": None,                      # [z]
            "CMvec": None,                  # [x,y,z] center of mass
        },
        "NRM_dir_stat": {
            "DANG": None,                   # = 5.587 (one number only)
            "NRMdev": None,                 # = 5.22  (one number only)
            "dir_vec_free": None            # [x, y, z]
        },
        "alpha_stat": {
            "alpha": None,                  # = 4.550 (one number only)
        },
        "theta_stat": {
            "theta": None,                  # = 38.66 (one number only)
        },
        "gamma_stat": {
            "gamma": None,                  # = 2.544 (one number only)
        }
    },

    "check_statistics": {
        "n_ptrm_check_stat": {
            "n_ptrm_check": None,           # = 3 (one integer only)
        },
        "max_and_cum_ptrm_check_stat": {
            "check_perc": None,             # = 10.01 (one number only)
            "d_CK": None,                   # = 5.65 (one number only)
            "DRAT": None,                   # = 4.98 (one number only)
            "maxDEV": None,                 # = 7.11 (one number only)
            "CDRAT": None,                  # = 5.01 (one number only)
            "CDRAT_prime": None,            # = 6.34 (one number only)
            "DRATS": None,                  # = 5.92 (one number only)
            "DRATS_prime": None,            # = 7.49 (one number only)
            "Mean_DRAT": None,              # = 1.67 (one number only)
            "Mean_DRAT_prime": None,        # = 2.11 (one number only)
            "Mean_DEV": None,               # = 2.38 (one number only)
            "Mean_DEV_prime": None,         # = 3.01 (one number only)
            "L": None,                      # = 7.50e-6 (one number only)
        },

        "dpal_ptrm_check_stat": {
            "d_pal": None,                  # = 10.01 (one number only)
        },

        "Tail_check_stat": {
            "n_tail_check": None,           # = 0 (one integer only)
            "DRAT_Tail": None,              # = x.xx (one number only)
            "dTR": None,                    # = x.xx (one number only)
            "MD_VDS": None,                 # = x.xx (one number only)
            "dt_star": None,                # = x.xx (one number only)
        },

        "Additivity_check": {
            "n_add": None,                  # = 0 (one integer only)
            "d_AC": None,                   # = x.xx (one number only)
        },
    },

    "anisotropy_statistics": {
        "Aniso_tensor": {
            "s_tensor": None,                # (one number only)
        },
        "Anisotropy_Correction": {
            "aniso_c": None,                # (one number only)
            "Banc_aniso_corr": None,        # (one number only)
        },
        "aniso_alteration": {
            "dTRM_anis": None,              # (one number only)
        },
    },

    "curv_arai_statistics": {
        "AraiCurvature": {
            "k_prime": None,               # (one number only)
            "SSE_k_prime": None,           # (one number only)
        }
    },

    "SCAT_parameter": {
        "SCAT_stat": {
            "a1": None,                     # (one number only)
            "a2": None,                     # (one number only)
            "s1": None,                     # (one number only)
            "s2": None,                     # (one number only)
            "SCAT": None,                   # = 1 ( or = 0 ) (one number only)
        }
    },

    "best_fit_lines": {
        "best_fit_line_Arai": {
            "l_start": None,                # [x, y]
            "l_end": None,                  # [x, y]
        },
        "SCAT_box_Arai": {
            "l1_start": None,               # [x, y]
            "l1_end": None,                 # [x, y]
            "l2_start": None,               # [x, y]
            "l2_end": None,                 # [x, y]
        },
        "best_fit_lines_Zijderveld": {
            "line_H_UpN": None,             # [[x1,y1], [x2,y2]] (Horizontal North/UP points for line in Zijderveld plot)
            "line_V_UpN": None,             # [[x1,y1], [x2,y2]] (Vertical North/UP points for line in Zijderveld plot)
            "line_H_UpW": None,             # [[x1,y1], [x2,y2]] (Horizontal West/UP points for line in Zijderveld plot)
            "line_V_UpW": None,             # [[x1,y1], [x2,y2]] (Vertical West/UP points for line in Zijderveld plot)
        }
    }
}

################################################################################################################################################################

MSP_suitcase = {
    "input": {
        "site": None,                       # data table with all specimens used
        "selection": None,                  # data table with selected selection for bootstrap
        "alpha": None,                      # (one number) user input between 0.0 - 1.0
    },

    "MSP_Q_calc": {
        "Q_DB": None,                       # [['name', x-coord, y-coord], [...], .. ]
        "Q_DSC": None,                      # [['name', x-coord, y-coord], [...], .. ]
        "Eps_alt": None,                    # [['name', number], [...], .. ] percent with 1.00 (100%) BELANGRIJK
    },

    "MSP_Q_calc_corr": {
        "Q_DB": None,                       # [['name', x-coord, y-coord], [...], .. ]
        "Q_DSC": None,                      # [['name', x-coord, y-coord], [...], .. ]
        "Eps_alt": None,                    # [['name', number], [...], .. ] percent with 1.00 (100%) BELANGRIJK
    },

    "MSP_results_Q_DB": {
        "PI": None,                         # (one number only)
        "avg_eps_alt": None,                # (one number only)
        "delta_b": None,                    # (one number only)
        "r_sq": None,                       # (one number only)
        "chi_sq": None,                     # (one number only)
        "Line_fig": None,                   # [[x1,y1], [x2,y2]] line through point 1 and 2
    },

    "MSP_results_Q_DSC": {
        "PI": None,                         # (one number only)
        "avg_eps_alt": None,                # (one number only)
        "delta_b": None,                    # (one number only)
        "r_sq": None,                       # (one number only)
        "chi_sq": None,                     # (one number only)
        "Line_fig": None,                   # [[x1,y1], [x2,y2]] line through point 1 and 2
    },

    "MSP_results_Q_DSC_corr": {
        "PI": None,                         # (one number only)
        "avg_eps_alt": None,                # (one number only)
        "delta_b": None,                    # (one number only)
        "r_sq": None,                       # (one number only)
        "chi_sq": None,                     # (one number only)
        "Line_fig": None,                   # [[x1,y1], [x2,y2]] line through point 1 and 2
    }
}

################################################################################################################################################################
pseudo_Thellier_suitcase = {
    "input": {
        "specimen": None,
        "specimen_selection": None
    },

    "preprocessed": {
        "msrmnts": {
            "zerofield_steps": None,
            "zerofield_steps_all": None,
            "ARM_acq_steps": None,
            "ARM_acq_steps_all": None,
            "ARM_dem_steps": None,
        },

        "basics": {
            "nrm0": None,                   # = 1.04e-5 (one number only)
            "ptrm": None,                   # [{step, x, y, z}]
            "NRM_rem": None,                # [{step, x, y, z}]
            "NRM_vec_all": None,            # [x,y,z]
            "ptrm_vec_all": None,           # [x,y,z]
            "y_nrm": None,                  # [y]
            "xy_temp": None,                # [step]
            "NRM_vec_select": None,         # [x,y,z]
            "x_ptrm": None,                 # [x]
            "x_temp": None,                 # [step]
            "ptrm_gained_vec": None,        # [x, y, z]
            "xBar": None,                   # = 3.04e-6 (one number only)
            "yBar": None,                   # = 1.04e-6 (one number only)
        },

        "basics_pTh": {
            "ARM_rem": None,                # [{step, x, y, z}]
            "ARMd_vec_all": None,           # [x,y,z]
            "y_arm_d": None,                # [y]
            "arm_step": None,                # [step]
            "ARMd_vec": None,               # [x,y,z]
            "yBar_arm": None,               # = 1.04e-6 (one number only)
        },

        "field_basics": {
            "Blab": None,                   # = 50 (one number only)
            "field_dir_vec_initial": None,  # = x, y, z (one line only)
            "field_dir_vec": None,          # = x, y, z (one line only)
        },
    },

    "arai_statistics": {
        "PI_est": {
            "n": None,                      # = 8 (one number only)
            "b_slope": None,                # = 1.02 (one number only)
            "b_sigma": None,                # = 0.04 (one number only)
            "b_beta": None,                 # = 1.02 (one number only)
            "SE_anc": None,                  # = 1.02 (one number only)
            "sumx": None,                  # = 1.02 (one number only)
            "sumy": None,                  # = 1.02 (one number only)
            "sum_xy": None,                  # = 1.02 (one number only)
        },
        "intercept_stats": {
            "Yint": None,                   # = 6.02e-6 (one number only)
            "Xint": None                    # = 7.02e-6 (one number only)
        },
        "VDS_stats": {
            "VDS": None,                    # = 1.12e-5 (one number only)
            "FRAC": None,                   # = 0.42 (one number only)
            "GAP_MAX": None                 # = 0.20 (one number only)
        },
        "xy_prime_stat": {
            "x_prime": None,                # [x]
            "y_prime": None,                # [y]
        },
        "delta_prime_stat": {
            "delta_x_prime": None,          # = 5.258e-6 (one number only)
            "delta_y_prime": None,          # = 5.352e-6 (one number only)
        },
        "f_stat": {
            "f": None,                      # = 0.795 (one number only)
        },
        "f_VDS_stat": {
            "f_VDS": None,                  # = 0.475 (one number only)
        },
        "g_stat": {
            "g": None,                      # = 0.857 (one number only)
            "g_lim": None,                  # = 0.818 (one number only)
        },
        "q_stat": {
            "q": None,                      # = 13.7 (one number only)
        },
        "w_stat": {
            "w": None,                      # = (one number only)
        },
        "rsq_stat": {
            "rsq_corr": None,               # = (one number only)
            "rsq_det": None,                # = (one number only)
        },
    },

    "directional_statistics": {
        "mean_dir_stat": {
            "Mdec_free": None,              # = 88.81  (one number only)
            "Minc_free": None,              # = -51.33 (one number only)
            "MAD_free": None,               # = 6.26   (one number only)
            "Mdec_anc": None,               # = 96.01  (one number only)
            "Minc_anc": None,               # = -50.829(one number only)
            "MAD_anc": None,                # = 3.483  (one number only)
            "x": None,                      # [x]
            "y": None,                      # [y]
            "z": None,                      # [z]
            "CMvec": None,                  # [x,y,z] center of mass
        },
        "NRM_dir_stat": {
            "DANG": None,                   # = 5.587 (one number only)
            "NRMdev": None,                 # = 5.22  (one number only)
            "dir_vec_free": None            # [x, y, z]
        },
        "alpha_stat": {
            "alpha": None,                  # = 4.550 (one number only)
        },
        "theta_stat": {
            "theta": None,                  # = 38.66 (one number only)
        },
        "gamma_stat": {
            "gamma": None,                  # = 2.544 (one number only)
        }
    },

    "curv_arai_statistics": {
        "AraiCurvature": {
            "k_prime": None,                # (one number only)
            "SSE_k_prime": None,            # (one number only)
        }
    },


    "pTh_statistics": {

        "initial": {
            "Yint_AA": None,                     # (one number only)
            "Yint_DD": None,                     # (one number only)
            "b_slope_AA": None,                 # (one number only)
            "b_slope_DD": None,                 # (one number only)
            "rsq_corr_AA": None,                # (one number only)
            "rsq_corr_DD": None,                # (one number only)
        },
        "curvature_AA": {
            "k_prime_AA": None,                # (one number only)
            "SSE_k_prime_AA": None,            # (one number only)
        },
        "curvature_DD": {
            "k_prime_DD": None,                # (one number only)
            "SSE_k_prime_DD": None,            # (one number only)
        },
    },

    "SCAT_parameter": {
        "SCAT_stat": {
            "a1": None,                     # (one number only)
            "a2": None,                     # (one number only)
            "s1": None,                     # (one number only)
            "s2": None,                     # (one number only)
            "SCAT": None,                   # = 1 ( or = 0 ) (one number only)
        }
    },

    "best_fit_lines": {
        "best_fit_line_Arai": {
            "l_start": None,                # [x, y]
            "l_end": None,                  # [x, y]
        },
        "best_fit_line_pTh_AA": {
            "l_start_AA": None,                # [x, y]
            "l_end_AA": None,                  # [x, y]
        },
        "best_fit_line_pTh_DD": {
            "l_start_DD": None,                # [x, y]
            "l_end_DD": None,                  # [x, y]
        },
        "SCAT_box_Arai": {
            "l1_start": None,               # [x, y]
            "l1_end": None,                 # [x, y]
            "l2_start": None,               # [x, y]
            "l2_end": None,                 # [x, y]
        },
        "best_fit_lines_Zijderveld": {
            "line_H_UpN": None,             # [[x1,y1], [x2,y2]] (Horizontal North/UP points for line in Zijderveld plot)
            "line_V_UpN": None,             # [[x1,y1], [x2,y2]] (Vertical North/UP points for line in Zijderveld plot)
            "line_H_UpW": None,             # [[x1,y1], [x2,y2]] (Horizontal West/UP points for line in Zijderveld plot)
            "line_V_UpW": None,             # [[x1,y1], [x2,y2]] (Vertical West/UP points for line in Zijderveld plot)
        },
    },
}
