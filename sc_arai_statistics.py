import math
import helpers


"""
This script contains the functions for the Arai-plot statistics, according to the SPD (Standard Paleointensity Definitions).

Paterson, G. A., L. Tauxe, A. J. Biggin, R. Shaar, and L. C. Jonestrask (2014), On improving the selection of Thellier-type paleointensity data, Geochem. Geophys. Geosyst., doi: 10.1002/2013GC005135
"""


def PI_est(sc):
    """
    Function to calculate the parameters to get the paleointensity estimate. The number of points used for the estimate of the best-fit line (n), the slope of the best-fit line on the Arai plot (b_slope), the relative scatter of the best-fit line (b_beta), standard error of the slope (b_sigma), and the standard error (SE_anc) on the paleointensity estimate.

    input: from the preprocessed data obtain the x and y values for Arai plot, also the xBar and yBar. from the field basics import the lab strength, Blab
    output: n, b_slope, b_sigma, b_beta, SE_anc, sumx, sumy, sum_xy
    """
    # input:    preprocessed/       basics[x_ptrm, y_nrm, xBar, yBar]
    #           preprocessed/       field_basics[Blab]
    # output:   arai_statistics/    PI_est[n, b_slope, b_sigma, b_beta, SE_anc, sumx, sumy, sum_xy]
    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    xBar = sc["preprocessed"]["basics"]["xBar"]
    yBar = sc["preprocessed"]["basics"]["yBar"]

    Blab = sc["preprocessed"]["field_basics"]["Blab"]

    n = len(x_ptrm)

    # Part (1) of b_slope equation
    sum_xy = 0

    for i in range(0, len(x_ptrm)):
        sum_xy += (x_ptrm[i] - xBar) * (y_nrm[i] - yBar)

    if sum_xy < 0:
        sign = -1
    elif sum_xy > 0:
        sign = 1
    else:
        sign = 0

    # part (2) of b_slope equation sumx en sumy
    sumx = 0
    sumy = 0
    for i in range(0, len(x_ptrm)):
        sumx += (x_ptrm[i] - xBar)**2
        sumy += (y_nrm[i] - yBar)**2

    # part(1) * part(2) gives b_slope
    b_slope = sign * math.sqrt(sumy / sumx)
    b_sigma = math.sqrt((2 * sumy - 2 * (b_slope) * sum_xy) / ((len(x_ptrm) - 2) * sumx))  # without absolute bars around the b_slope
    b_beta = b_sigma / abs(b_slope)
    SE_anc = b_sigma * Blab         # associated standard error

    sc["arai_statistics"]["PI_est"]["n"] = n
    sc["arai_statistics"]["PI_est"]["b_slope"] = b_slope
    sc["arai_statistics"]["PI_est"]["b_sigma"] = b_sigma
    sc["arai_statistics"]["PI_est"]["b_beta"] = b_beta
    sc["arai_statistics"]["PI_est"]["SE_anc"] = SE_anc
    sc["arai_statistics"]["PI_est"]["sumx"] = sumx
    sc["arai_statistics"]["PI_est"]["sumy"] = sumy
    sc["arai_statistics"]["PI_est"]["sum_xy"] = sum_xy
    return sc


def PI_Banc_est(sc):
    """
    Function to calculate the paleointensity estimate from the slope of the Arai plot and the strength of the labfield

    input: lab field strength (Blab), slope of the best-fit line (b_slope)
    output: paleointensity estimate (B_anc)
    """
    # input:    preprocessed/       field_basics[Blab]
    #           arai_statistics/    PI_est[b_slope]
    # output:   arai_statistics/    PI_Banc_est[B_anc]
    Blab = sc["preprocessed"]["field_basics"]["Blab"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]

    # from slope and sigma to paleointensity estimate
    B_anc = abs(b_slope) * Blab     # paleointensity estimate

    sc["arai_statistics"]["PI_Banc_est"]["B_anc"] = B_anc
    return sc


def intercept_stats(sc):
    """
    Function to calculate the intercept parameters. The y-axis, NRM, intercept of the best-fit line on the Arai plot (Yint). And the x-axis, TRM, intercept of the best-fit line on the Arai plot (Xint).

    input: slope of the best-fit line (b_slope), the mean TRM and NRM values of the selected data (xBar and yBar)
    output: intercepts of the best-fit line (Yint and Xint)
    """
    # input:    preprocessed/       basics[xBar, yBar, b_slope]
    #           arai_statistics/    PI_est[b_slope]
    # output:   arai_statistics/    intercept_stats[Yint, Xint]
    xBar = sc["preprocessed"]["basics"]["xBar"]
    yBar = sc["preprocessed"]["basics"]["yBar"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]

    Yint = 0
    Xint = 0
    Yint = yBar - b_slope * xBar
    Xint = -1 * (Yint) / b_slope

    sc["arai_statistics"]["intercept_stats"]["Xint"] = Xint
    sc["arai_statistics"]["intercept_stats"]["Yint"] = Yint
    return sc


def VDS_stats(sc):
    """
    Vector difference sum (VDS) calculations. These use the entire NRM vector. calculated here: the VDS, FRAC and GAP-MAX. FRAC and GAP-MAX proposed by (Shaar and Tauxe, 2013). FRAC is the NRM fraction used for the best-fit line. GAP-MAX is the maximum gap between two points.

    input: The entire NRM vector and the NRM vector of the selection
    output: VDS, FRAC, GAP_MAX
    """
    # input:    preprocessed/       basics[NRM_vec_all, NRM_vec_select]
    # output:   arai_statistics/    VDS_stats[VDS, FRAC, GAP_MAX]
    NRM_vec_select = sc["preprocessed"]["basics"]["NRM_vec_select"]
    NRM_vec_all = sc["preprocessed"]["basics"]["NRM_vec_all"]

    # calculate difference vector, use a helper function

    def calc_diff(MagVec):
        diff = []
        for i in range(len(MagVec) - 1):
            diff.append(helpers.difference(MagVec[i], MagVec[i + 1]))
        return diff

    diff_total_NRM = calc_diff(NRM_vec_all)  # this is only the second part of the equation, need to add NRMmax

    # calculate norm of the last NRM step, the NRM_max
    NRM_max = helpers.norm(NRM_vec_all[len(NRM_vec_all) - 1])

    diff_total_NRM.append(NRM_max)  # append NRM max tot the VDS vector
    VDS = sum(diff_total_NRM)  # take the sum of all elements

    # Calculation of FRAC
    # calculate the difference vector for the correct selection
    # NRM_vec_select is the vector only containing the NRM of the selection of interest
    diff_select_NRM = calc_diff(NRM_vec_select)

    FRAC = 0
    FRAC = sum(diff_select_NRM) / VDS

    # Calculation of GAP-MAX
    GAP_MAX = 0
    GAP_MAX = max(diff_select_NRM) / sum(diff_select_NRM)

    sc["arai_statistics"]["VDS_stats"]["VDS"] = VDS
    sc["arai_statistics"]["VDS_stats"]["FRAC"] = FRAC
    sc["arai_statistics"]["VDS_stats"]["GAP_MAX"] = GAP_MAX
    return sc


def xy_prime_stat(sc):
    """
    Function to calculate the x' and y' parameters. The primes are the x and y points on the Arai plot projected on to the best-fit line.

    input: all x and y points (x_ptrm & y_nrm), the slope of the best-fit line (b_slope), and the y-intercept (Yint)
    output: x' and y' (x_prime, y_prime)
    """
    # input:    preprocessed/       basics[x_ptrm, y_nrm]
    #           arai_statistics/    PI_est[b_slope]
    #                               intercept_stats[Yint]
    # output:   arai_statistics/    xy_prime_stat[x_prime, y_prime]

    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]

    x_prime = []
    y_prime = []

    for i in range(0, len(x_ptrm)):
        x_prime.append(0.5 * (x_ptrm[i] + (y_nrm[i] - Yint) / b_slope))
        y_prime.append(0.5 * (y_nrm[i] + b_slope * x_ptrm[i] + Yint))

    sc["arai_statistics"]["xy_prime_stat"]["x_prime"] = x_prime
    sc["arai_statistics"]["xy_prime_stat"]["y_prime"] = y_prime
    return sc


def delta_prime_stat(sc):
    """
    Function to calculate the TRM and NRM lengths of the best-fit line on the Arai plot, delta_x_prime and delta_y_prime.

    input: x' and y' (x_prime, y_prime)
    output: delta_x_prime and delta_y_prime
    """
    # input:    arai_statistics/    xy_prime_stat[x_prime, y_prime]
    # output:   arai_statistics/    delta_prime_stat[delta_x_prime, delta_y_prime]
    x_prime = sc["arai_statistics"]["xy_prime_stat"]["x_prime"]
    y_prime = sc["arai_statistics"]["xy_prime_stat"]["y_prime"]

    delta_x_prime = 0
    delta_y_prime = 0

    delta_x_prime = abs(max(x_prime) - min(x_prime))
    delta_y_prime = abs(max(y_prime) - min(y_prime))

    sc["arai_statistics"]["delta_prime_stat"]["delta_x_prime"] = delta_x_prime
    sc["arai_statistics"]["delta_prime_stat"]["delta_y_prime"] = delta_y_prime
    return sc


def f_stat(sc):
    """
    Function to calculate the NRM fraction used for the best-fit on an Arai diagram (Coe et al., 1978).

    input: the y-intercept (Yint) and the NRM length of the best-fit line on the Arai plot (delta_y_prime)
    output: the NRM fraction (f)
    """
    # input:    arai_statistics/    intercept_stats[Yint]
    #           arai_statistics/    delta_prime_stat[delta_y_prime]
    # output:   arai_statistics/    f_stat[f]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]
    delta_y_prime = sc["arai_statistics"]["delta_prime_stat"]["delta_y_prime"]

    f = delta_y_prime / abs(Yint)

    sc["arai_statistics"]["f_stat"]["f"] = f
    return sc


def f_VDS_stat(sc):
    """
    Function to calculate the NRM fraction used for the best-fit on an Arai diagram calculated as a vector difference sum (Tauxe and Staudigel, 2004)

    input: Vector difference sum (VDS), NRM length of the best-fit line on the Arai plot (delta_y_prime)
    output: NRM fraction (f_VDS)
    """
    # input:    arai_statistics/    VDS_stats[VDS]
    #           arai_statistics/    delta_prime_stat[delta_y_prime]
    # output:   arai_statistics/    f_VDS_stat[f_VDS]
    VDS = sc["arai_statistics"]["VDS_stats"]["VDS"]
    delta_y_prime = sc["arai_statistics"]["delta_prime_stat"]["delta_y_prime"]

    f_VDS = delta_y_prime / VDS

    sc["arai_statistics"]["f_VDS_stat"]["f_VDS"] = f_VDS
    return sc


def g_stat(sc):
    """
    Function to calculate the gap factor (g), and the upper limit of g (g_lim). The gap reflects the average spacing of the selected Arai plot points along the best-fit line.

    input: number of points (n), y points on the Arai plot projected on to the best-fit line (y_prime), NRM length of the best-fit line on the Arai plot (delta_y_prime)
    output: (g) and (g_lim)
    """
    # input:    arai_statistics/    PI_est[n]
    #           arai_statistics/    xy_prime_stat[y_prime]
    #           arai_statistics/    delta_prime_stat[delta_y_prime]
    # output:   arai_statistics/    g_stat[g, g_lim]
    n = sc["arai_statistics"]["PI_est"]["n"]
    y_prime = sc["arai_statistics"]["xy_prime_stat"]["y_prime"]
    delta_y_prime = sc["arai_statistics"]["delta_prime_stat"]["delta_y_prime"]

    g_sq_diff = []
    for i in range(len(y_prime) - 1):
        g_sq_diff.append((y_prime[i + 1] - y_prime[i])**2)

    g = 1.0 - (sum(g_sq_diff) / (delta_y_prime**2))
    g_lim = (n - 2.0) / (n - 1.0)

    sc["arai_statistics"]["g_stat"]["g"] = g
    sc["arai_statistics"]["g_stat"]["g_lim"] = g_lim
    return sc


def q_stat(sc):
    """
    Function to calculate the quality factor (q), (Coe et al., 1978). q is a measure of the overall quality of the paleointensity estimate.

    input: the slope of the best-fit line, (b_slope), relative scatter of the best-fit line (b_beta), standard error of the slope (b_sigma), the NRM fraction (f), and the gap factor (g).
    output: the quality factor (q)
    """
    # input:    arai_statistics/    PI_est[b_slope, b_beta, b_sigma]
    #           arai_statistics/    f_stat[f]
    #           arai_statistics/    g_stat[g]
    # output:   arai_statistics/    q_stat[q]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]
    b_beta = sc["arai_statistics"]["PI_est"]["b_beta"]
    b_sigma = sc["arai_statistics"]["PI_est"]["b_sigma"]
    f = sc["arai_statistics"]["f_stat"]["f"]
    g = sc["arai_statistics"]["g_stat"]["g"]

    q = (abs(b_slope) * f * g) / b_sigma

    # also possible with b_beta, should give the same answer
    q_beta = (f * g) / b_beta

    sc["arai_statistics"]["q_stat"]["q"] = q
    return sc


def w_stat(sc):
    """
    Function to calculate the weighting factor (w) of Prévot et al. (1985).

    input: quality factor (q), the number of points used for the estimate of the best-fit line (n)
    output: weighting factor (w)
    """
    # input:    arai_statistics/    q_stat[q]
    #           arai_statistics/    PI_est[n]
    # output:   arai_statistics/    w_stat[w]
    q = sc["arai_statistics"]["q_stat"]["q"]
    n = sc["arai_statistics"]["PI_est"]["n"]

    w = q / math.sqrt(n - 2)

    sc["arai_statistics"]["w_stat"]["w"] = w
    return sc


def rsq_stat(sc):
    """
    Function to calculate the correlation coefficient (rsq_corr), and the coefficient of determination (rsq_det)

    input: sumx, sumy, sum_xy, y_prime, y_nrm
    output: rsq_corr, rsq_det
    """
    # input:    arai_statistics/    PI_est[sumx, sumy, sum_xy]
    #                               xy_prime_stat[y_prime]
    #           preprocessed/       basics[y_nrm]
    # output:   arai_statistics/    rsq_stat[rsq_corr, rsq_det]
    sumx = sc["arai_statistics"]["PI_est"]["sumx"]
    sumy = sc["arai_statistics"]["PI_est"]["sumy"]
    sum_xy = sc["arai_statistics"]["PI_est"]["sum_xy"]

    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    y_prime = sc["arai_statistics"]["xy_prime_stat"]["y_prime"]

    rsq_corr = (sum_xy)**2 / (sumx * sumy)

    sumyy_prime = 0
    for i in range(0, len(y_nrm)):
        sumyy_prime += (y_nrm[i] - y_prime[i])**2

    rsq_det = 1 - (sumyy_prime / sumy)

    sc["arai_statistics"]["rsq_stat"]["rsq_corr"] = rsq_corr
    sc["arai_statistics"]["rsq_stat"]["rsq_det"] = rsq_det
    return sc
