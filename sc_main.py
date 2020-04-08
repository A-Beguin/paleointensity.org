import sc_config
import helpers
import helpers_curv

import sc_preprocessing
import sc_arai_statistics
import sc_directional
import sc_anisotropy
import sc_checks
import sc_bestfit_lines
import sc_curvature_arai
import sc_SCAT_param

import sc_pTh_statistics

import sc_MSP_results

import math

"""
This is the main body to create and fill the suitcases. All functions are imported and are called in the correct order in the pipeline. There are separate functions for thermal_thellier, microwave_thellier, pseudo_thellier, and MSP
"""


def thermal_thellier(specimen, selection):
    # create and fill the suitcase
    suitcase = sc_config.new_Thellier_suitcase()
    suitcase["input"]["specimen"] = specimen
    suitcase["input"]["specimen_selection"] = selection

    # define the pipeline
    pipeline = [
        sc_preprocessing.splitup_msrmnts,
        sc_preprocessing.prep_anisotropy_tensor,
        sc_preprocessing.calc_nrm_ptrm,
        sc_preprocessing.basics,
        sc_preprocessing.field_basics,
        sc_preprocessing.checks,

        sc_arai_statistics.PI_est,
        sc_arai_statistics.PI_Banc_est,
        sc_arai_statistics.intercept_stats,
        sc_arai_statistics.VDS_stats,
        sc_arai_statistics.xy_prime_stat,
        sc_arai_statistics.delta_prime_stat,
        sc_arai_statistics.f_stat,
        sc_arai_statistics.f_VDS_stat,
        sc_arai_statistics.g_stat,
        sc_arai_statistics.q_stat,
        sc_arai_statistics.w_stat,
        sc_arai_statistics.rsq_stat,

        sc_directional.mean_dir_stat,
        sc_directional.NRM_dir_stat,
        sc_directional.alpha_stat,
        sc_directional.theta_stat,
        sc_directional.gamma_stat,

        sc_anisotropy.s_tensor_calc,
        sc_anisotropy.anisotropy_calc,
        sc_anisotropy.aniso_alteration_calc,

        sc_checks.n_ptrm_check_stat,
        sc_checks.max_and_cum_ptrm_check_stat,
        sc_checks.dpal_ptrm_check_stat,
        sc_checks.Tail_check_stat,
        sc_checks.Tail_check_dtstar_stat,
        sc_checks.Additivity_check_stat,

        # sc_curvature_arai_old.AraiCurvature,
        sc_curvature_arai.AraiCurvature,

        sc_SCAT_param.SCAT_stat,

        sc_bestfit_lines.best_fit_line_Arai,
        sc_bestfit_lines.SCAT_box_Arai,
        sc_bestfit_lines.best_fit_lines_Zijderveld,

    ]

    suitcase = exec_pipeline(pipeline, suitcase)

    return suitcase


def microwave_thellier(specimen, selection):
    return thermal_thellier(specimen, selection)


def pseudo_thellier(specimen, selection):
    # create and fill the suitcase
    suitcase = sc_config.new_pseudo_Thellier_suitcase()
    suitcase["input"]["specimen"] = specimen
    suitcase["input"]["specimen_selection"] = selection

    # define the pipeline
    pipeline = [
        sc_preprocessing.splitup_msrmnts,
        sc_preprocessing.calc_nrm_ptrm,
        sc_preprocessing.basics,
        sc_preprocessing.basics_pTh,
        sc_preprocessing.field_basics_pTh,

        # # sc_preprocessing.checks,                     # NOT necessary for pTH

        sc_arai_statistics.PI_est,
        sc_arai_statistics.intercept_stats,
        sc_arai_statistics.VDS_stats,
        sc_arai_statistics.xy_prime_stat,
        sc_arai_statistics.delta_prime_stat,
        sc_arai_statistics.f_stat,
        sc_arai_statistics.f_VDS_stat,
        sc_arai_statistics.g_stat,
        sc_arai_statistics.q_stat,
        sc_arai_statistics.w_stat,
        sc_arai_statistics.rsq_stat,

        sc_directional.mean_dir_stat,
        sc_directional.NRM_dir_stat,
        sc_directional.alpha_stat,
        sc_directional.theta_stat,
        sc_directional.gamma_stat,

        sc_curvature_arai.AraiCurvature,

        sc_pTh_statistics.rsq_corr_pTh_stat,
        sc_pTh_statistics.k_prime_DD_pTh_stat,
        sc_pTh_statistics.k_prime_AA_pTh_stat,

        sc_SCAT_param.SCAT_stat,

        sc_bestfit_lines.best_fit_line_Arai,
        sc_bestfit_lines.best_fit_line_pTh_AA,
        sc_bestfit_lines.best_fit_line_pTh_DD,
        sc_bestfit_lines.SCAT_box_Arai,
        sc_bestfit_lines.best_fit_lines_Zijderveld,
    ]

    suitcase = exec_pipeline(pipeline, suitcase)

    return suitcase


def MSP(site, selection, alpha, NumCycles, cutOff):
    # create and fill the suitcase
    suitcase = sc_config.new_MSP_suitcase()
    suitcase["input"]["site"] = site
    suitcase["input"]["selection"] = selection
    suitcase["input"]["alpha"] = alpha
    suitcase["input"]["NumCycles"] = NumCycles
    suitcase["input"]["cutOff"] = cutOff

    # define the pipeline
    pipeline = [
        sc_MSP_results.MSP_Q_calc,
        sc_MSP_results.MSP_results_calc,
    ]

    suitcase = exec_pipeline(pipeline, suitcase)

    return suitcase


def exec_pipeline(pipeline, suitcase):
    for step in pipeline:
        suitcase = step(suitcase)

    return suitcase
