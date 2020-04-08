import math
import helpers
import helpers_curv

"""
This script contains the functions for the pTh statistics. This is for the additional two plots, the demag-demag and ARM-ARM plot.
"""


def rsq_corr_pTh_stat(sc):
    """
    Function to calculate the additional pTh statistics for the two plots, demag-demag and ARM-ARM.

    input: arm_step, y_arm_d, ARMd_vec, yBar_arm, xBar, yBar, sumx, sumy, x_ptrm, y_nrm, xy_temp, x_temp, sumx, sumy
    output: Yint_AA, Yint_DD, b_slope_AA, b_slope_DD, rsq_corr_AA, rsq_corr_DD, f_resid

    """

    # input:    preprocessed/       msrmnts[ARM_dem_steps]
    #                               basics_pTh[arm_step, y_arm_d, ARMd_vec, yBar_arm]
    #                               basics[xBar, yBar, sumx, sumy, x_ptrm, y_nrm, xy_temp, x_temp]
    #           arai_statistics/    PI_est[sumx, sumy]
    # output:   pTh_statistics/     initial[Yint_AA, Yint_DD, b_slope_AA, b_slope_DD, rsq_corr_AA, rsq_corr_DD, f_resid]

    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]

    # only do all calculations if the ARM_dem is given, this implies that the ARM_aq is also given
    if len(ARM_dem_steps) != 0:

        arm_step = sc["preprocessed"]["basics_pTh"]["arm_step"]
        y_arm_d = sc["preprocessed"]["basics_pTh"]["y_arm_d"]
        ARMd_vec = sc["preprocessed"]["basics_pTh"]["ARMd_vec"]
        yBar_arm = sc["preprocessed"]["basics_pTh"]["yBar_arm"]

        xBar = sc["preprocessed"]["basics"]["xBar"]
        yBar = sc["preprocessed"]["basics"]["yBar"]
        sumx = sc["arai_statistics"]["PI_est"]["sumx"]
        sumy = sc["arai_statistics"]["PI_est"]["sumy"]
        x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
        y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
        xy_temp = sc["preprocessed"]["basics"]["xy_temp"]
        x_temp = sc["preprocessed"]["basics"]["x_temp"]

        # calculate the sumy for the arm demag
        sumy_arm = 0
        for i in range(len(y_arm_d)):
            sumy_arm += (y_arm_d[i] - yBar_arm)**2

        # first for ARM-ARM
        # on x axis is now x_ptrm with xBar -> sumx
        # on y axis is the y_arm_d with yBar_arm -> sumy_arm

        sum_xy_AA = 0
        # check is the steps are the same and calculate sum_xy_AA
        for i in range(len(y_arm_d)):
            for j in range(len(x_temp)):
                if arm_step[i] == x_temp[j]:
                    sum_xy_AA += (x_ptrm[i] - xBar) * (y_arm_d[i] - yBar_arm)

        rsq_corr_AA = (sum_xy_AA)**2 / (sumx * sumy_arm)

        # repeat for Demag-Demag
        # on x axis is now y_nrm with yBar -> sumy
        # on y axis is again the y_arm_d with yBar_arm -> sumy_arm

        sum_xy_DD = 0
        # check is the steps are the same and calculate sum_xy_AA
        for i in range(len(y_arm_d)):
            for j in range(len(xy_temp)):
                if arm_step[i] == xy_temp[j]:
                    sum_xy_DD += (y_nrm[i] - yBar) * (y_arm_d[i] - yBar_arm)

        rsq_corr_DD = (sum_xy_DD)**2 / (sumy * sumy_arm)

        # calculate the slope for the ARM-ARM plot

        if sum_xy_AA < 0:
            sign = -1
        elif sum_xy_AA > 0:
            sign = 1
        else:
            sign = 0

        # part (2) of b_slope equation sumx en sumy_arm is already calculated

        # part(1) * part(2) gives b_slope
        b_slope_AA = sign * math.sqrt(sumy_arm / sumx)

        # caclulate the y intercept for AA plot
        Yint_AA = yBar_arm - b_slope_AA * xBar

        # calculate the slope for the demag-demag plot

        if sum_xy_DD < 0:
            sign = -1
        elif sum_xy_DD > 0:
            sign = 1
        else:
            sign = 0

        # part (2) of b_slope equation sumx en sumy_arm is already calcluated

        # part(1) * part(2) gives b_slope
        b_slope_DD = sign * math.sqrt(sumy_arm / sumy)

        # calculate the y intercept for DD plot
        Yint_DD = yBar_arm - b_slope_DD * yBar

        y_prime = []
        for i in range(0, len(y_arm_d)):
            y_prime.append(0.5 * (y_arm_d[i] + b_slope_DD * y_nrm[i] + Yint_DD))

        delta_y_prime = abs(max(y_prime) - min(y_prime))

        f_resid = abs(Yint_DD) / delta_y_prime

        sc["pTh_statistics"]["initial"]["Yint_AA"] = Yint_AA
        sc["pTh_statistics"]["initial"]["Yint_DD"] = Yint_DD
        sc["pTh_statistics"]["initial"]["b_slope_AA"] = b_slope_AA
        sc["pTh_statistics"]["initial"]["b_slope_DD"] = b_slope_DD
        sc["pTh_statistics"]["initial"]["rsq_corr_AA"] = rsq_corr_AA
        sc["pTh_statistics"]["initial"]["rsq_corr_DD"] = rsq_corr_DD
        sc["pTh_statistics"]["initial"]["f_resid"] = f_resid

    return sc


def k_prime_DD_pTh_stat(sc):
    """
    Function to calculate the curvature for the demag-demag plot.

    input: the x and y points for the demag-demag plot, that is the y_nrm and y_arm_d
    output: curvature k and SEE
    """
    # input:    preprocessed/   msrmnts[ARM_dem_steps]
    #                           basics[y_nrm]
    #                           basics_pTh[y_arm_d]
    # output:   pTh_statistics/ curvature_DD/   k_prime_DD[k]
    #                                           SSE_k_prime_DD[SSE]

    # only do all calculations if the ARM_dem is given, this implies that the ARM_aq is also given
    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]

    if len(ARM_dem_steps) != 0:
        x = sc["preprocessed"]["basics"]["y_nrm"]
        y = sc["preprocessed"]["basics_pTh"]["y_arm_d"]

        if len(x) > 3:
            (k, SSE) = helpers_curv.kprime_calc(x, y)

            sc["pTh_statistics"]["curvature_DD"]["k_prime_DD"] = k
            sc["pTh_statistics"]["curvature_DD"]["SSE_k_prime_DD"] = SSE
    return(sc)


def k_prime_AA_pTh_stat(sc):
    """
    Function to calculate the curvature for the ARM-ARM plot.

    input: the x and y points for the ARM-ARM plot, that is the x_ptrm and y_arm_d
    output: curvature k and SEE
    """
    # input:    preprocessed/   msrmnts[ARM_dem_steps]
    #                           basics[x_ptrm]
    #                           basics_pTh[y_arm_d]
    # output:   pTh_statistics/ curvature_AA/   k_prime_AA[k]
    #                                           SSE_k_prime_AA[SSE]

    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]

    if len(ARM_dem_steps) != 0:
        x = sc["preprocessed"]["basics"]["x_ptrm"]
        y = sc["preprocessed"]["basics_pTh"]["y_arm_d"]

        if len(x) > 3:
            (k, SSE) = helpers_curv.kprime_calc(x, y)

            sc["pTh_statistics"]["curvature_AA"]["k_prime_AA"] = k
            sc["pTh_statistics"]["curvature_AA"]["SSE_k_prime_AA"] = SSE
    return(sc)
