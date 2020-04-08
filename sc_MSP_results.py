import math
import helpers


"""
This script contains two functions for the calculation of the MSP fractions (MSP_Q_calc) and results in the suitcase (MSP_results_calc). The calculations follow the calculations from the MSP-Tool by Monster et al., 2015.

Monster, M. W. L., de Groot, L. V, & Dekkers, M. J. (2015). MSP-Tool: A VBA-Based Software Tool for the Analysis of Multispecimen Paleointensity Data. Frontiers in Earth Science, 3, 3–9. https://doi.org/10.3389/feart.2015.00086
"""


def MSP_Q_calc(sc):
    """
    Function to calculate the MSP fractions (Q). Using the input selection and the alpha as defined by the input.

    input: selection of data (selection), alpha
    output: Q_DB, Q_DSC, Eps_alt and the corrected versions
    """

    # input:    input[selection, alpha]
    # output:   MSP_Q_calc[Q_DB, Q_DSC, Eps_alt]
    #           MSP_Q_calc_corr[Q_DB, Q_DSC, Eps_alt]

    selection = sc["input"]["selection"]
    alpha = sc["input"]["alpha"]

    # split in measurements m0 m1 m2 m3 m4 with multiple specimens per list
    m0 = list(filter(lambda m: m['type'] == 0, selection))
    m1 = list(filter(lambda m: m['type'] == 1, selection))
    m2 = list(filter(lambda m: m['type'] == 2, selection))
    m3 = list(filter(lambda m: m['type'] == 3, selection))
    m4 = list(filter(lambda m: m['type'] == 4, selection))

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
                Eps_alt.append([name, None])

            else:  # calculate also the Q_DSC ratio and all the parameters
                Q_DB.append([name, H_lab, (m_m1 - m_m0) / m_m0])  # no corrected version, only "normal version"

                # first check for the corrected version of un-corrected version for Q_DSC & parameter calculations
                if (c_str == "_corr"):
                    # calculate corrections
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
                    m_m2 = helpers.norm(NRMrem) - helpers.norm(m2pTRM)  # exception to the rule
                    m_m3 = helpers.norm(NRMrem) + helpers.norm(m3pTRM)
                    m_m4 = helpers.norm(NRMrem) + helpers.norm(m4pTRM)

                Q_DSC.append([name, H_lab, 2 * ((1 + alpha) * m_m1 - m_m0 - alpha * m_m3) / (2 * m_m0 - m_m1 - m_m2)])

                Eps = (m_m4 - m_m1) / m_m1
                Eps_alt.append([name, Eps])

        sc["MSP_Q_calc" + c_str]["Q_DB"] = Q_DB
        sc["MSP_Q_calc" + c_str]["Q_DSC"] = Q_DSC
        sc["MSP_Q_calc" + c_str]["Eps_alt"] = Eps_alt

    params_noCorr_corr("", alpha, m0, m1, m2, m3, m4)
    params_noCorr_corr("_corr", alpha, m0, m1, m2, m3, m4)

    return sc


def MSP_results_calc(sc):
    """
    Function to calculate the MSP results for the MSP fraction (Q) for the different options (DB, DSC, DSC_corr) using a function "results" with the input for results being the corresponding Q factor.

    input: MSP fractions Q_DB, Q_DSC, Q_DSC_corr
    output: for the different options (DB, DSC, DSC_corr) the paleointensity estimate (PI), average epsilon alteration (avg_eps_alt), the intercept difference (delta_b), r-squared (r_sq), chi-squared (chi_sq), (x,y) coordinated for two point of the best-fit line (Line_fig)
    """
    # input:    MSP_Q_calc[Q_DB, Q_DSC, Eps_alt]
    #           MSP_Q_calc_corr[Q_DSC, Eps_alt]
    # output:   MSP_results_Q_DB[PI, avg_eps_alt, delta_b, r_sq, chi_sq]
    #           MSP_results_Q_DSC[PI, avg_eps_alt, delta_b, r_sq, chi_sq]
    #           MSP_results_Q_DSC_corr[PI, avg_eps_alt, delta_b, r_sq, chi_sq]

    # get the results for both Q_DB & Q_DSC with new function & for the corrected values
    def results(Q_str, c_str, Eps_alt):
        Q = sc["MSP_Q_calc" + c_str]["Q_" + Q_str]

        x = []
        y = []
        EpsAlt = []

        # if (Q_str == "DSC"):  # check for None in the specimen list
        for i in range(len(Q)):
            if Q[i][1] != None:
                x.append(Q[i][1])
                y.append(Q[i][2])
                EpsAlt.append(Eps_alt[i][1])

        N = len(x)  # sumber of specimens

        if N > 1:  # if N> 1 then you have enough specimens to calculate the Linear regression
            # determine x and y coordinates, and Eps

            # calculate Sx, Sy, Sxx, Syy, Sxy
            Sx = sum(x)
            Sy = sum(y)

            Sxy = helpers.dot_product(x, y)
            Sxx = helpers.dot_product(x, x)
            Syy = helpers.dot_product(y, y)

            # calculate linear regression coefficients
            LRb = (N * Sxy - Sx * Sy) / (N * Sxx - Sx**2)
            LRa = Sy / N - LRb * Sx / N

            # determine PI
            PI = -1 * LRa / LRb

            # get two points for the linear regression line
            x1 = -1
            x2 = 500
            y1 = LRa + LRb * x1
            y2 = LRa + LRb * x2
            Line_fig = [[x1, y1], [x2, y2]]

            # calculate the average, agv X and Y for r squared

            avg_x = Sx / N
            avg_y = Sy / N

            xDiff = helpers.list_min_num(x, avg_x)
            yDiff = helpers.list_min_num(y, avg_y)

            x2Sum = helpers.dot_product(xDiff, xDiff)
            y2Sum = helpers.dot_product(yDiff, yDiff)
            xySum = helpers.dot_product(xDiff, yDiff)

            r_sq = (xySum / math.sqrt(x2Sum * y2Sum))**2

            # difficlt expression: yexp = Hlab[i]*LRb - LRa
            Yexp = helpers.list_min_num(helpers.list_mult_num(x, LRb), -1 * LRa)

            yminYexp = helpers.list_min_list(y, Yexp)

            ChiSum = helpers.dot_product(yminYexp, yminYexp)
            chi_sq = ChiSum / N

            # calculate the average epsilon for DSC only

            if (Q_str == "DB"):
                delta_b = None
                avg_eps_alt = None
            else:
                delta_b = LRa + 1
                avg_eps_alt = sum(EpsAlt) / N

            sc["MSP_results_Q_" + Q_str + c_str]["PI"] = PI
            sc["MSP_results_Q_" + Q_str + c_str]["avg_eps_alt"] = avg_eps_alt
            sc["MSP_results_Q_" + Q_str + c_str]["delta_b"] = delta_b
            sc["MSP_results_Q_" + Q_str + c_str]["r_sq"] = r_sq
            sc["MSP_results_Q_" + Q_str + c_str]["chi_sq"] = chi_sq
            sc["MSP_results_Q_" + Q_str + c_str]["Line_fig"] = Line_fig  # line through point 1 and 2, [[x1,y1], [x2,y2]]

    # find the correct epsilon
    Eps_alt = sc["MSP_Q_calc"]["Eps_alt"]
    Eps_alt_corr = sc["MSP_Q_calc_corr"]["Eps_alt"]

    # run the different options
    results("DB", "", Eps_alt)
    results("DSC", "", Eps_alt)
    results("DSC", "_corr", Eps_alt_corr)

    return sc
