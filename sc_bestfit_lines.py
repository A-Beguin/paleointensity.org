import math
import helpers

"""
This script contains the functions for the calculation of the best-fit lines for the Arai-plot, SCAT-boxes, Zijdervel-plot, and pTh-plots.
"""


def best_fit_line_Arai(sc):
    """
    Function to calculate the best-fit line in the Arai-plot, this uses the pTRM for the x-axis and the y-intercept for the y-axis.

    input: x_ptrm, nrm0, b_slope, Yint
    output: (x,y) coordinates for the start (l_start) and end (l_end) of the best fit line
    """
    # input:    preprocessed/       basics[x_ptrm, nrm0]
    #           arai_statistics/    PI_est[b_slope]
    #           arai_statistics/    intercept_stats[Yint]
    # output:   best_fit_lines/     best_fit_line_Arai[l_start, l_end]

    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    nrm0 = sc["preprocessed"]["basics"]["nrm0"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]

    l_start = []  # x and y coordinate of start best fit line
    l_end = []  # x and y coordinate of end best fit line

    # x coordinates of best fit line start and end at the first and last ptrm step
    lx0 = x_ptrm[0]
    lx1 = x_ptrm[len(x_ptrm) - 1]

    # calculate corresponding y coordinate for the best fit line
    ly0 = b_slope * lx0 + Yint
    ly1 = b_slope * lx1 + Yint

    # normalize by nrm0
    l_start = [lx0 / nrm0, ly0 / nrm0]
    l_end = [lx1 / nrm0, ly1 / nrm0]

    sc["best_fit_lines"]["best_fit_line_Arai"]["l_start"] = l_start
    sc["best_fit_lines"]["best_fit_line_Arai"]["l_end"] = l_end
    return sc


def SCAT_box_Arai(sc):
    """
    Function to calculate the upper and lower SCAT lines, using the pre-calculated boundingbox points.

    input: information on the x-axis (x_ptrm) and the NRM0 (nrm0) for plotting. The SCAT boundingbox points (a1,a2,s1,s2)
    output: the (x,y) coordinates for the start and end of line 1 (l1_start, l1_end), and for line 2 (l2_start, l2_end)
    """
    # input:    preprocessed/   basics[x_ptrm, nrm0]
    #           SCAT_parameter/ SCAT_stat[a1,a2,s1,s2]
    # output:   best_fit_lines/ SCAT_box_Arai[l1_start, l1_end, l2_start, l2_end]
    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    nrm0 = sc["preprocessed"]["basics"]["nrm0"]
    a1 = sc["SCAT_parameter"]["SCAT_stat"]["a1"]
    a2 = sc["SCAT_parameter"]["SCAT_stat"]["a2"]
    s1 = sc["SCAT_parameter"]["SCAT_stat"]["s1"]
    s2 = sc["SCAT_parameter"]["SCAT_stat"]["s2"]

    l1_start = []  # x and y coordinate of start line 1
    l1_end = []  # x and y coordinate of end line 1

    l2_start = []  # x and y coordinate of start line 2
    l2_end = []  # x and y coordinate of end line 2

    # x coordinates of two SCAT_Box lines start and end at the first and last ptrm step
    x_start = x_ptrm[0]
    x_end = x_ptrm[-1]

    # The y coordinates are calculated with the two line formulas and the x start and x end coordinates
    # line_1 = a1 - (a1/s2)x        &      line_2 = a2 - (a2/s1)x
    l1_y_start = a1 - (a1 / s2) * x_start
    l1_y_end = a1 - (a1 / s2) * x_end

    l2_y_start = a2 - (a2 / s1) * x_start
    l2_y_end = a2 - (a2 / s1) * x_end

    l1_start = [x_start / nrm0, l1_y_start / nrm0]
    l1_end = [x_end / nrm0, l1_y_end / nrm0]

    l2_start = [x_start / nrm0, l2_y_start / nrm0]
    l2_end = [x_end / nrm0, l2_y_end / nrm0]

    sc["best_fit_lines"]["SCAT_box_Arai"]["l1_start"] = l1_start
    sc["best_fit_lines"]["SCAT_box_Arai"]["l1_end"] = l1_end
    sc["best_fit_lines"]["SCAT_box_Arai"]["l2_start"] = l2_start
    sc["best_fit_lines"]["SCAT_box_Arai"]["l2_end"] = l2_end
    return sc


def best_fit_lines_Zijderveld(sc):
    """
    Function to calculate the direction of the Best - Fit lines for the Zijderveld, use the Free floating direction

    input: Free-floating direction (Mdec_free, Minc_free) and the center-of-mass vector (CMvec) to calculate the size of the plot
    output: The (x,y) coordinates for the horizontal (line_H_) and vertical (line_V_) lines for the Zijderveld diagram, for the two options, Up-North (UpN) or Up-West (UpW)
    """
    # input:    directional_statistics/ mean_dir_stat[CMvec, Mdec_free, Minc_free ]
    # output:   best_fit_lines/         best_fit_lines_Zijderveld[line_H_UpN, line_V_UpN, line_H_UpW, line_V_UpW]

    CMvec = sc["directional_statistics"]["mean_dir_stat"]["CMvec"]
    Mdec_free = sc["directional_statistics"]["mean_dir_stat"]["Mdec_free"]
    Minc_free = sc["directional_statistics"]["mean_dir_stat"]["Minc_free"]

    # first find in which order of magnitude the point for the lines should be, they should be outside the area of the zijderveld plot
    # find the maximum absolute value for the center of mass and take an order of magnitude bigger to be sure to outside of the plot
    Abs_CMvec = []
    for i in range(len(CMvec)):
        Abs_CMvec.append(abs(CMvec[i]))

    M = max(Abs_CMvec) * 10

    # get the direction vector with the correct magnitude M
    Dvec = helpers.dir2cart(Mdec_free, Minc_free, M)

    # Center of mass - Direction vector is the first point
    P1 = helpers.list_min_list(CMvec, Dvec)

    # Center of mass + Direction vector is the second point
    P2 = helpers.list_plus_list(CMvec, Dvec)

    # North Up
    line_H_UpN = [[P1[1], P1[0]], [P2[1], P2[0]]]
    line_V_UpN = [[P1[1], -1 * P1[2]], [P2[1], -1 * P2[2]]]

    # West Up
    line_H_UpW = [[P1[0], -1 * P1[1]], [P2[0], -1 * P2[1]]]
    line_V_UpW = [[P1[0], -1 * P1[2]], [P2[0], -1 * P2[2]]]

    sc["best_fit_lines"]["best_fit_lines_Zijderveld"]["line_H_UpN"] = line_H_UpN
    sc["best_fit_lines"]["best_fit_lines_Zijderveld"]["line_V_UpN"] = line_V_UpN

    sc["best_fit_lines"]["best_fit_lines_Zijderveld"]["line_H_UpW"] = line_H_UpW
    sc["best_fit_lines"]["best_fit_lines_Zijderveld"]["line_V_UpW"] = line_V_UpW
    return sc


def best_fit_line_pTh_AA(sc):
    """
    Function to calculate the best-fit line for the pseudo-Thellier ARM-ARM plot. Using the ARM steps (x_ptrm) and the ARM-Demag steps (y_arm_d).

    input: ARM info (x_ptrm), ARM-Dem info (y_arm_d), the y-intercept for the ARM-ARM plot (Yint_AA), and the slope of the ARM-ARM plot (b_slope_AA)
    output: (x,y) coordinates for the start (l_start) and end (l_end) of the best-fit line
    """
    # input:    preprocessed/       basics[x_ptrm, nrm0]
    #                               basics_pTh[y_arm_d]
    #           pTh_statistics/     initial[Yint_AA, b_slope_AA]
    # output:   best_fit_lines/     best_fit_line_pTh_AA[l_start, l_end]

    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]

    if len(ARM_dem_steps) != 0:
        x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
        y_arm_d = sc["preprocessed"]["basics_pTh"]["y_arm_d"]
        nrm0 = sc["preprocessed"]["basics"]["nrm0"]
        b_slope = sc["pTh_statistics"]["initial"]["b_slope_AA"]
        Yint = sc["pTh_statistics"]["initial"]["Yint_AA"]

        l_start = []  # x and y coordinate of start best fit line
        l_end = []  # x and y coordinate of end best fit line

        # x coordinates of best fit line start and end at the first and last ptrm step
        lx0 = x_ptrm[0]
        lx1 = x_ptrm[len(x_ptrm) - 1]

        # calculate corresponding y coordinate for the best fit line
        ly0 = b_slope * lx0 + Yint
        ly1 = b_slope * lx1 + Yint

        # normalize by nrm0
        l_start = [lx0 / nrm0, ly0 / nrm0]
        l_end = [lx1 / nrm0, ly1 / nrm0]

        sc["best_fit_lines"]["best_fit_line_pTh_AA"]["l_start_AA"] = l_start
        sc["best_fit_lines"]["best_fit_line_pTh_AA"]["l_end_AA"] = l_end
    return sc


def best_fit_line_pTh_DD(sc):
    """
    Function to calculate the best-fit line for the pseudo-Thellier demag-demag plot. Using the NRM steps (y_nrm) and the ARM-Demag steps (y_arm_d).

    input: NRM info (y_nrm, nrm0), ARM-Dem info (y_arm_d), the y-intercept for the demag-demag plot (Yint_DD), and the slope of the demag-demag plot (b_slope_DD)
    output: (x,y) coordinates for the start (l_start) and end (l_end) of the best-fit line
    """
    # input:    preprocessed/       basics[y_nrm, nrm0]
    #                               basics_pTh[y_arm_d]
    #           pTh_statistics/     initial[Yint_DD, b_slope_DD]
    # output:   best_fit_lines/     best_fit_line_pTh_DD[l_start, l_end]

    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]

    if len(ARM_dem_steps) != 0:
        y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
        y_arm_d = sc["preprocessed"]["basics_pTh"]["y_arm_d"]
        nrm0 = sc["preprocessed"]["basics"]["nrm0"]
        b_slope = sc["pTh_statistics"]["initial"]["b_slope_DD"]
        Yint = sc["pTh_statistics"]["initial"]["Yint_DD"]

        l_start = []  # x and y coordinate of start best fit line
        l_end = []  # x and y coordinate of end best fit line

        # x coordinates of best fit line start and end at the first and last ptrm step
        lx0 = y_nrm[0]
        lx1 = y_nrm[len(y_nrm) - 1]

        # calculate corresponding y coordinate for the best fit line
        ly0 = b_slope * lx0 + Yint
        ly1 = b_slope * lx1 + Yint

        # normalize by nrm0
        l_start = [lx0 / nrm0, ly0 / nrm0]
        l_end = [lx1 / nrm0, ly1 / nrm0]

        sc["best_fit_lines"]["best_fit_line_pTh_DD"]["l_start_DD"] = l_start
        sc["best_fit_lines"]["best_fit_line_pTh_DD"]["l_end_DD"] = l_end
    return sc
