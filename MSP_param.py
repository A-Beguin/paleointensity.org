import math
import helpers

"""

This script is used to calculate the MSP parameters. This is independent of the suitcase.

"""


def MSP_params_calc(data, alpha):
    """
    This function is used to get the MSP parameters, and contains one functions: params_noCorr_corr . This function returns the parameters; Q_DB, Q_DSC, mu_ds, H_max, H_est, Eps_alt, Eps_alt_abs, Err_alt, Err_ds, Err_total. The function is called twice, once for the alignment corrected parameters and once for non-corrected parameters.

    input: data for the specimens and alpha used for DSC calculations.
    output: params & params_corr
    """
    # input:    input           [data, alpha]
    #
    # output:   MSP_params[Q_DB, Q_DSC, mu_ds, H_max, H_est, Eps_alt, Eps_alt_abs, Err_alt, Err_ds, Err_total]
    # output:   MSP_params_corr[Q_DB, Q_DSC, mu_ds, H_max, H_est, Eps_alt, Eps_alt_abs, Err_alt, Err_ds, Err_total]

    m0 = list(filter(lambda m: m['type'] == 0, data))
    m1 = list(filter(lambda m: m['type'] == 1, data))
    m2 = list(filter(lambda m: m['type'] == 2, data))
    m3 = list(filter(lambda m: m['type'] == 3, data))
    m4 = list(filter(lambda m: m['type'] == 4, data))

    def params_noCorr_corr(c_str, alpha, m0, m1, m2, m3, m4):
        Q_DB = []
        Q_DSC = []
        mu_ds = []
        H_max = []
        H_est = []
        Eps_alt = []
        Eps_alt_abs = []
        Err_alt = []
        Err_ds = []
        Err_total = []
        Err_alt_abs = []

        num_specimens = len(m0)
        for i in range(num_specimens):

            name = m0[i]["specimen"]  # specimen name is the same for m0-m4
            H_lab = m1[i]["lab_field"]  # the x-axis coordinate is the AF field of the m1-m4 steps

            # do the calculations
            m_m0 = m0[i]["total_m"]
            m_m1 = m1[i]["total_m"]
            m_m2 = m2[i]["total_m"]
            m_m3 = m3[i]["total_m"]
            m_m4 = m4[i]["total_m"]

            if (m_m2 == None) or (m_m3 == None) or (m_m4 == None):
                Q_DB.append([name, H_lab, (m_m1 - m_m0) / m_m0])
                Q_DSC.append([name, None, None])
                mu_ds.append([name, None])
                H_max.append([name, None])
                H_est.append([name, None])
                Eps_alt.append([name, None])
                Eps_alt_abs.append([name, None])
                Err_alt.append([name, None])
                Err_ds.append([name, None])
                Err_total.append([name, None])

            else:  # calculate also the Q_DSC ratio and all the parameters

                Q_DB.append([name, H_lab, (m_m1 - m_m0) / m_m0])  # no corrected version, only "normal version"

                #  check for the corrected version of un-corrected version for Q_DSC & prameter calculations
                if (c_str == "_corr"):
                    m0M = [m0[i]["x"], m0[i]["y"], m0[i]["z"]]
                    m1M = [m1[i]["x"], m1[i]["y"], m1[i]["z"]]
                    m2M = [m2[i]["x"], m2[i]["y"], m2[i]["z"]]
                    m3M = [m3[i]["x"], m3[i]["y"], m3[i]["z"]]
                    m4M = [m4[i]["x"], m4[i]["y"], m4[i]["z"]]

                    NRMrem = helpers.list_mult_num(helpers.list_plus_list(m1M, m2M), 0.5)

                    m1pTRM = helpers.list_min_list(m1M, NRMrem)
                    m2pTRM = helpers.list_min_list(m2M, NRMrem)
                    m3pTRM = helpers.list_min_list(m3M, NRMrem)
                    m4pTRM = helpers.list_min_list(m4M, NRMrem)

                    m_m0 = m0[i]["total_m"]                                 # m_m0_corr
                    m_m1 = helpers.norm(NRMrem) + helpers.norm(m1pTRM)      # m_m1_corr
                    m_m2 = helpers.norm(NRMrem) - helpers.norm(m2pTRM)      # exception to the rule
                    m_m3 = helpers.norm(NRMrem) + helpers.norm(m3pTRM)
                    m_m4 = helpers.norm(NRMrem) + helpers.norm(m4pTRM)

                Q_DSC.append([name, H_lab, 2 * ((1 + alpha) * m_m1 - m_m0 - alpha * m_m3) / (2 * m_m0 - m_m1 - m_m2)])
                mu_ds.append([name, (m_m1 - m_m3) / (m_m3 - 0.5 * (m_m1 + m_m2))])
                H_max.append([name, (2 * m_m0 - m_m1 - m_m2) / (m_m1 - m_m2) * H_lab])
                H_est.append([name, (2 * m_m0 - m_m1 - m_m2) / ((1 + 2 * alpha) * m_m1 - 2 * alpha * m_m3 - m_m2) * H_lab])

                Eps = (m_m4 - m_m1) / m_m1

                Eps_alt.append([name, Eps])
                Eps_alt_abs.append([name, abs((m_m1 - m_m4) / m_m1)])

                # calculate the error estimates
                # nummerator & denominator
                num = 2 * ((1 + alpha) * m_m1 - m_m0 - alpha * m_m3)
                den = (2 * m_m0) - m_m1 - m_m2

                # partial derivatives of Q_DSC
                d_num_m1 = 2 * (1 + alpha)
                d_num_m2 = 0
                d_num_m3 = -2 * alpha

                d_den_m1 = -1
                d_den_m2 = -1
                d_den_m3 = 0

                # terms
                Term_1 = m_m1 * ((den * d_num_m1) - (num * d_den_m1)) / (den)**2
                Term_2 = m_m2 * ((den * d_num_m2) - (num * d_den_m2)) / (den)**2
                Term_3 = m_m3 * ((den * d_num_m3) - (num * d_den_m3)) / (den)**2

                dQ_DSC_alt = Eps**2 * (Term_1**2 + Term_2**2 + Term_3**2)

                dQ_DSC_ds = (((m_m3 - m_m1) / den)**2) / 3

                Err_alt.append([name, dQ_DSC_alt])
                Err_ds.append([name, dQ_DSC_ds])
                Err_total.append([name, dQ_DSC_alt + dQ_DSC_ds])

        return [Q_DB, Q_DSC, mu_ds, H_max, H_est, Eps_alt, Eps_alt_abs, Err_alt, Err_ds, Err_total]

    params = params_noCorr_corr("", alpha, m0, m1, m2, m3, m4)
    params_corr = params_noCorr_corr("_corr", alpha, m0, m1, m2, m3, m4)

    return [params, params_corr]
