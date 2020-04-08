import math

"""
This script contains the functions for the SCAT statistic of the Arai-plot, according to the SPD (Standard Paleointensity Definitions).

Paterson, G. A., L. Tauxe, A. J. Biggin, R. Shaar, and L. C. Jonestrask (2014), On improving the selection of Thellier-type paleointensity data, Geochem. Geophys. Geosyst., doi: 10.1002/2013GC005135
"""


def SCAT_stat(sc):
    """
    Function to calculate the SCAT statistic, proposed by Shaar and Tauxe (2013).

    input: all points that are within the selection of the best-fit line. This also includes pTRM-checks, tail-checks
    output: SCAT and the four coordinates for the scat-lines
    """

    # input:    preprocessed/       basics[x_ptrm, y_nrm, xBar, yBar]
    #                               msrmnts[ptrm_check_steps,tail_check_steps  ]
    #                               checks[x_ptrm_check, y_ptrm_check, y_tail_check, x_tail_check, x_temp_ptrm_check, SCAT_ptrm_check_step, y_temp_tail_check, SCAT_tail_check_step]
    #           arai_statistics/    PI_est[b_slope]
    # output:   SCAT_parameter/     SCAT_stat[SCAT, a1, a2, s1, s2]
    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    xBar = sc["preprocessed"]["basics"]["xBar"]
    yBar = sc["preprocessed"]["basics"]["yBar"]

    check_ptrm = sc["preprocessed"]["msrmnts"]["ptrm_check_steps"]
    check_tail = sc["preprocessed"]["msrmnts"]["tail_check_steps"]

    # only try to read in the x & y of the ptrm and the tail check is these are defined in the mrmnts
    # this is added to the script so that the pseudo Thellier code can use the same SCAT function
    if len(check_ptrm) > 0:
        x_ptrm_check = sc["preprocessed"]["checks"]["x_ptrm_check"]
        y_ptrm_check = sc["preprocessed"]["checks"]["y_ptrm_check"]
        x_temp_ptrm_check = sc["preprocessed"]["checks"]["x_temp_ptrm_check"]
        SCAT_ptrm_check_step = sc["preprocessed"]["checks"]["SCAT_ptrm_check_step"]
    if len(check_tail) > 0:
        y_tail_check = sc["preprocessed"]["checks"]["y_tail_check"]
        x_tail_check = sc["preprocessed"]["checks"]["x_tail_check"]
        y_temp_tail_check = sc["preprocessed"]["checks"]["y_temp_tail_check"]
        SCAT_tail_check_step = sc["preprocessed"]["checks"]["SCAT_tail_check_step"]

    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]

    all_point_x = []
    all_point_y = []

    # SCAT is determined for all point that are within the selection. the pTRM checks previous to the lowest step for the best fit should not be taken into account. These are within the x_ptrm_check steps but not in y_ptrm_check. extra parameter to check this is SCAT_ptrm_check_step
    # same holds for tail checks

    # find all the x and corresponding y coordinates for all the points in Arai plot
    # first append all the NRM and pTRM gained points
    for l in range(0, len(x_ptrm)):
        all_point_x.append(x_ptrm[l])
        all_point_y.append(y_nrm[l])

    # if ptrm checks and tail checks exist, also append these
    if len(check_ptrm) > 0:
        # second append x and y coordinates of the pTRM checks
        for l in range(0, len(SCAT_ptrm_check_step)):
            # the x_ptrm_check and y_ptrm_check are not necessarily the same length. If an ptrm check is performed before the selection of the best fit line the len(x_ptrm_check) > len(y_ptrm_check) for SCAT only take the ptrm checks that fall within the best fit selection. Use SCAT_ptrm_check_step & find corresponding steps for both x an y
            if len(x_ptrm_check) == len(y_ptrm_check):
                all_point_x.append(x_ptrm_check[l])
                all_point_y.append(y_ptrm_check[l])
            else:
                all_point_y.append(y_ptrm_check[l])
                for i in range(len(x_ptrm_check)):
                    if x_temp_ptrm_check[i] == SCAT_ptrm_check_step[l]:
                        all_point_x.append(x_ptrm_check[i])

    if len(check_tail) > 0:  # x_tail_check != None:
        # third append x and y coordinates of the tail-checks
        for l in range(0, len(SCAT_tail_check_step)):
            if len(y_tail_check) == len(x_tail_check):
                all_point_x.append(x_tail_check[l])
                all_point_y.append(y_tail_check[l])
            else:
                all_point_x.append(x_tail_check[l])
                for i in range(len(y_tail_check)):
                    if y_temp_tail_check[i] == SCAT_tail_check_step[l]:
                        all_point_y.append(y_tail_check[i])

    num_point = (len(all_point_x))
    num_point_y = (len(all_point_y))     # should be the same as num_point, since equal number of x and y coordinates

    # do SCAT statistics
    # two lines:    y1= a1 - (abs(b_slope) + 2*sigma_th)*xBar
    #               y2= a2 - (abs(b_slope) - 2*sigma_th)*xBar
    beta_th = 0.1
    sigma_th = abs(b_slope) * beta_th

    # determine y-axis intercepts a1, a2
    a1 = yBar + (abs(b_slope) + 2 * sigma_th) * xBar
    a2 = yBar + (abs(b_slope) - 2 * sigma_th) * xBar

    # and determine x-axis intercepts s1, s2
    s1 = a1 / (abs(b_slope) + 2 * sigma_th)
    s2 = a2 / (abs(b_slope) - 2 * sigma_th)

    # difine two lines with the above axis-intercepts
    # line_1 = a1 - (a1/s2)x
    # line_2 = a2 - (a2/s1)x

    # for a point with (Xp,Yp) check SCAT, do this for all the points in Arai plot
    SCAT_count = 0
    for i in range(0, num_point):   # loop over each point in Arai plot and determine is it is within the SCAT Box
        Xp = all_point_x[i]
        Yp = all_point_y[i]

        l1_p = a1 - (a1 / s2) * Xp
        l2_p = a2 - (a2 / s1) * Xp

        if (Xp > s2):           # then SCAT should fail, so appoint 0
            SCAT_count += 0
        elif (Yp > l1_p):       # then SCAT should fail, so appoint 0
            SCAT_count += 0
        elif (Yp < l2_p):
            SCAT_count += 0     # then SCAT should fail, so appoint 0
        else:
            SCAT_count += 1     # else, SCAT is true for this point

    if (SCAT_count < num_point):    # if all points lie within the scat box, the SCAT_count should be equal to the num_point
        SCAT = 0                # SCAT fail
    else:
        SCAT = 1                # SCAT pass

    sc["SCAT_parameter"]["SCAT_stat"]["a1"] = a1
    sc["SCAT_parameter"]["SCAT_stat"]["a2"] = a2
    sc["SCAT_parameter"]["SCAT_stat"]["s1"] = s1
    sc["SCAT_parameter"]["SCAT_stat"]["s2"] = s2
    sc["SCAT_parameter"]["SCAT_stat"]["SCAT"] = SCAT
    return sc
