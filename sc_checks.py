import math
import helpers


"""
This script contains six functions to calculate pTRM-check, tail-check and additivity-check statistics. All according to the SPD, Paterson et al., 2014
"""

def n_ptrm_check_stat(sc):
    """
    Function that gives the number of pTRM checks used to analyze the best-fit segment on the Arai plot

    input: list of ptrm-check measurements
    output: n_ptrm_check
    """
    # input:    preprocessed/       checks[x_ptrm_check]
    # output:   check_statistics/   n_ptrm_check_stat[n_ptrm_check]
    x_ptrm_check = sc["preprocessed"]["checks"]["x_ptrm_check"]

    n_ptrm_check = len(x_ptrm_check)

    sc["check_statistics"]["n_ptrm_check_stat"]["n_ptrm_check"] = n_ptrm_check
    return sc


def max_and_cum_ptrm_check_stat(sc):
    """
    Function to calculate the maximum pTRM check parameters (5.2 SPD), and the cumulative pTRM check parameters (5.3 SPD).

    input: ptrm-check data (x_ptrm_check, x_temp_ptrm_check), information on the TRM data points (x_ptrm, x_temp, x_ptrm_all, x_temp_all), the x-axis intercept of the best-fit line (Xint), and delta_x_prime delta_y_prime
    output: check_perc, d_CK, DRAT, maxDEV, CDRAT, CDRAT_prime, DRATS, DRATS_prime, Mean_DRAT, Mean_DRAT_prime, Mean_DEV, Mean_DEV_prime, L
    """
    # input:    preprocessed/       checks[x_ptrm_check, x_temp_ptrm_check]
    #                               basics[x_ptrm, x_temp, x_ptrm_all, x_temp_all]
    #           arai_statistics/    intercept_stats[Xint]
    #                               delta_prime_stat[delta_x_prime, delta_y_prime]
    # output:   check_statistics/   max_and_cum_ptrm_check_stat[check_perc, d_CK, DRAT, maxDEV, CDRAT, CDRAT_prime, DRATS, DRATS_prime, Mean_DRAT, Mean_DRAT_prime, Mean_DEV, Mean_DEV_prime, L]
    x_ptrm_check = sc["preprocessed"]["checks"]["x_ptrm_check"]
    x_temp_ptrm_check = sc["preprocessed"]["checks"]["x_temp_ptrm_check"]
    x_ptrm = sc["preprocessed"]["basics"]["x_ptrm"]
    x_temp = sc["preprocessed"]["basics"]["x_temp"]

    Xint = sc["arai_statistics"]["intercept_stats"]["Xint"]
    delta_x_prime = sc["arai_statistics"]["delta_prime_stat"]["delta_x_prime"]
    delta_y_prime = sc["arai_statistics"]["delta_prime_stat"]["delta_y_prime"]

    x_ptrm_all = sc["preprocessed"]["basics"]["x_ptrm_all"]
    x_temp_all = sc["preprocessed"]["basics"]["x_temp_all"]

    check_perc_vec = []
    check_perc = []
    abs_delta_pTRM = []

    d_CK = []
    DRAT = []
    maxDEV = []

    delta_pTRM = []
    vec = []

    CDRAT = []
    CDRAT_prime = []
    DRATS = []
    DRATS_prime = []
    Mean_DRAT = []
    Mean_DRAT_prime = []
    Mean_DEV = []
    Mean_DEV_prime = []

    L = math.sqrt((delta_x_prime)**2 + (delta_y_prime)**2)

    n_ptrm_check = len(x_ptrm_check)
    if n_ptrm_check != 0:

        # statistic check(%) -> check_perc
        for i in range(len(x_ptrm_check)):
            for j in range(len(x_ptrm_all)):
                if x_temp_ptrm_check[i] == x_temp_all[j]:                       # For each check temp find the corresponding TRM temp
                    delta_pTRM.append(x_ptrm_check[i] - x_ptrm_all[j])          # scalar intensity diff of pTRMcheck with original pTRM
                    vec.append(j)                                           # vector with the index for pTRM vector, easier for rest of calculations

        for i in range(len(x_ptrm_check)):                                  # make a vector with the absolute diff produced by a pTRM check
            abs_delta_pTRM.append(abs(delta_pTRM[i]))
            check_perc_vec.append(abs_delta_pTRM[i] / x_ptrm_all[vec[i]] * 100) # normalized by TRM acquired at the corresponding heating step
        check_perc = max(check_perc_vec)                                    # find maximum of absolute difference vector

        # statistic d_CK
        max_delta_pTRM = max(abs_delta_pTRM)
        d_CK = (max_delta_pTRM / abs(Xint)) * 100

        # statistic DRAT
        DRAT = (max_delta_pTRM / L) * 100

        # statistic maxDEV
        maxDEV = (max_delta_pTRM / delta_x_prime) * 100

        # Cumulative parameters
        abs_sum_dpTRM = abs(sum(delta_pTRM))
        sum_abs_dpTRM = sum(abs_delta_pTRM)

        # statistic CDRAT & CDRAT_prime
        CDRAT = (abs_sum_dpTRM / L) * 100
        CDRAT_prime = (sum_abs_dpTRM / L) * 100

        # statistic DRATS & DRATS_prime
        DRATS = (abs_sum_dpTRM / x_ptrm[len(x_ptrm)-1]) * 100
        DRATS_prime = (sum_abs_dpTRM / x_ptrm[len(x_ptrm)-1]) * 100

        # statistic Mean_DRAT & Mean_DRAT_prime
        Mean_DRAT = (1. / len(x_ptrm_check)) * (abs_sum_dpTRM / L) * 100
        Mean_DRAT_prime = (1. / len(x_ptrm_check)) * (sum_abs_dpTRM / L) * 100

        # statistic Mean_DEV & Mean_DEV_prime
        Mean_DEV = (1. / len(x_ptrm_check)) * (abs_sum_dpTRM / delta_x_prime) * 100
        Mean_DEV_prime = (1. / len(x_ptrm_check)) * (sum_abs_dpTRM / delta_x_prime) * 100

    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["check_perc"] = check_perc
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["d_CK"] = d_CK
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["DRAT"] = DRAT
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["maxDEV"] = maxDEV
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["CDRAT"] = CDRAT
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["CDRAT_prime"] = CDRAT_prime
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["DRATS"] = DRATS
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["DRATS_prime"] = DRATS_prime
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["Mean_DRAT"] = Mean_DRAT
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["Mean_DRAT_prime"] = Mean_DRAT_prime
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["Mean_DEV"] = Mean_DEV
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["Mean_DEV_prime"] = Mean_DEV_prime
    sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["L"] = L
    return sc


def dpal_ptrm_check_stat(sc):
    """
    Function to calculate dpal statistic. A measure of cumulative alteration determined by the difference of the alteration corrected intensity estimate (Valet et al., 1996) and the uncorrected estimate, normalized by the uncorrected estimate (Leonhardt et al., 2004a).

    input: slope of the best-fit line, ptrm_checks, ptrm, y_nrm, yBar
    output: d_pal
    """
    # input:    preprocessed/       checks[ptrm_check]
    #           preprocessed/       basics[ptrm, y_nrm, yBar]
    #           arai_statistics/    PI_est[b_slope]
    # output:   check_statistics/   dpal_ptrm_check_stat[d_pal]

    ptrm_check = sc["preprocessed"]["checks"]["ptrm_check"]
    ptrm = sc["preprocessed"]["basics"]["ptrm"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    yBar = sc["preprocessed"]["basics"]["yBar"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]

    d_pal = []
    dT = []

    # first check is you have ptrm_checks in the data, if not do nothing
    if len(ptrm_check) != 0:
        # first determine d_pTRM, which is the vector difference between ptrm - ptrmcheck for each step, if no ptrm check performed this is 0, 0, 0
        s = 0
        sprev = 0
        for T in ptrm:       # look trough the calculated ptrm gained specimen data
            for C in ptrm_check:
                if (T["step"] == C["step"]):
                    dT.append([T["x"] - C["x"], T["y"] - C["y"], T["z"] - C["z"]])
                    s += 1
            if s <= sprev:
                dT.append([ 0, 0, 0])
            sprev = s

        # make the cumulative sum vector C
        C = []
        Ci = [0,0,0]
        for i in range(len(ptrm)):
            Ci = helpers.list_plus_list(Ci, dT[i])
            C.append(Ci)

        ptrm_list = []
        # add C to TRM
        for T in ptrm:
            ptrm_list.append([T["x"],T["y"],T["z"]])

        TRM_star = []
        for j in range(len(ptrm_list)):
            TRM_star.append(helpers.list_plus_list(ptrm_list[j],C[j]))

        x_ptrm_star = []
        for i in range(len(TRM_star)):
            x_ptrm_star.append(helpers.norm(TRM_star[i]))

        # calculate the "new slope"
        # copy form sc_arai_statiscics
        n = len(x_ptrm_star)
        xBar_star = sum(x_ptrm_star) / len(x_ptrm_star)

        # Part (1) of b_slope equation
        sum_xy = 0
        for i in range(0, len(x_ptrm_star)):
            sum_xy += (x_ptrm_star[i] - xBar_star) * (y_nrm[i] - yBar)

        if sum_xy < 0:
            sign = -1
        elif sum_xy > 0:
            sign = 1
        else:
            sign = 0

        # part (2) of b_slope equation sumx en sumy
        sumx = 0
        sumy = 0
        for i in range(0, len(x_ptrm_star)):
            sumx += (x_ptrm_star[i] - xBar_star)**2
            sumy += (y_nrm[i] - yBar)**2

        # part(1) * part(2) gives b_slope
        b_slope_star = sign * math.sqrt(sumy / sumx)

        ## stop copy
        d_pal = abs((b_slope - b_slope_star) / b_slope) * 100


    sc["check_statistics"]["dpal_ptrm_check_stat"]["d_pal"] = d_pal
    return sc


def Tail_check_stat(sc):
    """
    Function to calculate the pTRM tail checks statistic, excluding dt*.

    input: tail check info (y_tail_check, y_temp_tail_check), NRM info (y_nrm, xy_temp), other calculated statistics (VDS, L)
    output: n_tail_check, DRAT_Tail, dTR, MD_VDS
    """
    # input:    preprocessed/       checks[y_tail_check, y_temp_tail_check]
    #           preprocessed/       basics[y_nrm, xy_temp]
    #           arai_statistics/    intercept_stats[Yint]
    #           arai_statistics/    VDS_stats[VDS]
    #           check_statistics/   max_and_cum_ptrm_check_stat[L]
    # output:   check_statistics/   Tail_check_stat[n_tail_check, DRAT_Tail, dTR, MD_VDS]
    y_tail_check = sc["preprocessed"]["checks"]["y_tail_check"]
    y_temp_tail_check = sc["preprocessed"]["checks"]["y_temp_tail_check"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    xy_temp = sc["preprocessed"]["basics"]["xy_temp"]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]
    VDS = sc["arai_statistics"]["VDS_stats"]["VDS"]
    L = sc["check_statistics"]["max_and_cum_ptrm_check_stat"]["L"]

    y_nrm_all = sc["preprocessed"]["basics"]["y_nrm_all"]
    y_temp_all = sc["preprocessed"]["basics"]["y_temp_all"]

    DRAT_Tail = []
    dTR = []
    MD_VDS = []

    n_tail_check = len(y_tail_check)

    # y_tail_check can contain data from outside the selection data, use y_nmr_all and y_temp_all
    if n_tail_check != 0:
        delta_tail = []
        for i in range(len(y_tail_check)):
            for j in range(len(y_nrm_all)):
                if y_temp_tail_check[i] == y_temp_all[j]:              # For each check temp find the corresponding NRM temp
                    delta_tail.append(y_tail_check[i] - y_nrm_all[j])   # scalar intensity diff of tailCheck with original NRMrem

        abs_delta_tail = []
        for i in range(len(y_tail_check)):                          # make a vector with the absolute diff produced by a tail check
            abs_delta_tail.append(abs(delta_tail[i]))
        max_abs_dTail = max(abs_delta_tail)

        # statistic DRAT_Tail
        DRAT_Tail = (max_abs_dTail / L) * 100

        # statistic dTR
        dTR = (max_abs_dTail / abs(Yint)) * 100

        # statistic MD_VDS
        MD_VDS = (max_abs_dTail / VDS) * 100

    sc["check_statistics"]["Tail_check_stat"]["n_tail_check"] = n_tail_check
    sc["check_statistics"]["Tail_check_stat"]["DRAT_Tail"] = DRAT_Tail
    sc["check_statistics"]["Tail_check_stat"]["dTR"] = dTR
    sc["check_statistics"]["Tail_check_stat"]["MD_VDS"] = MD_VDS
    return sc


def Tail_check_dtstar_stat(sc):
    """
    Function to calculate the extent of a pTRM tail after correction for angular dependence (Leonhardt et al., 2004a; 2004b).

    input: y_tail_check, y_temp_tail_check, y_nrm, y_nrm_all, x_temp_all, y_temp_all, NRM_vec_select, NRM_vec_all, xy_temp, field_dir_vec, Yint, Xint, b_slope
    output: dt_star
    """
    # input:    preprocessed/       checks[y_tail_check, y_temp_tail_check]
    #                               basics[y_nrm, y_nrm_all, x_temp_all, y_temp_all, NRM_vec_select, NRM_vec_all, xy_temp]
    #                               field_basics[field_dir_vec]
    #           arai_statistics/    intercept_stats[Yint, Xint]
    #                               PI_est[b_slope]
    # output:   check_statistics/   Tail_check_stat[dt_star]
    tail_check_vec = sc["preprocessed"]["checks"]["tail_check_vec"]
    y_temp_tail_check = sc["preprocessed"]["checks"]["y_temp_tail_check"]
    y_nrm = sc["preprocessed"]["basics"]["y_nrm"]
    y_nrm_all = sc["preprocessed"]["basics"]["y_nrm_all"]
    x_temp_all = sc["preprocessed"]["basics"]["x_temp_all"]
    y_temp_all = sc["preprocessed"]["basics"]["y_temp_all"]
    xy_temp = sc["preprocessed"]["basics"]["xy_temp"]
    NRM_vec_select = sc["preprocessed"]["basics"]["NRM_vec_select"]
    NRM_vec_all = sc["preprocessed"]["basics"]["NRM_vec_all"]
    field_dir_vec = sc["preprocessed"]["field_basics"]["field_dir_vec"]
    b_slope = sc["arai_statistics"]["PI_est"]["b_slope"]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]
    Xint = sc["arai_statistics"]["intercept_stats"]["Xint"]

    dt_star = []
    t_star = []


    if len(y_temp_tail_check)!= 0:  # then you have at least one tailcheck, do calculations
        # make an index matrix for the steps that step tail = step nrm_all
        ind = []
        for i in range(len(y_temp_tail_check)):             # number of tail checks
            for j in range(len(y_temp_all)):                # loop over all NRM steps
                if y_temp_tail_check[i] == y_temp_all[j]:   # if step of tail check equals step op NRM, write index to "ind" this is used below
                    ind.append(j)

        # do calculations per tail check step
        for i in range(len(y_temp_tail_check)):
            MDx = tail_check_vec[i][0]
            MDy = tail_check_vec[i][1]
            MDz = tail_check_vec[i][2]

            Nx = NRM_vec_all[ind[i]][0]            # take the NRM vector with same step as tail_check
            Ny = NRM_vec_all[ind[i]][1]
            Nz = NRM_vec_all[ind[i]][2]


            a = helpers.list_div_num(NRM_vec_all[ind[i]],y_nrm_all[ind[i]]) # divide the NRM list by the corresponding moment
            theta_dt = helpers.get_angle_diff(a, field_dir_vec)             # get angle between nomalized NRM vector and field direction

            # calculate dH and dZ and inc_diff
            # according to orientation of Blab
            if abs(field_dir_vec[0]/helpers.norm(field_dir_vec)) ==1: # then allong x
                dH = math.sqrt(Ny**2 + Nz**2) - math.sqrt(MDy**2 + MDz**2)
                dZ = Nx - MDx
                Dec_F, Inc_F, R_F = helpers.cart2dir(field_dir_vec[0],field_dir_vec[1],field_dir_vec[2])
                Dec_N, Inc_N, R_N = helpers.cart2dir(Nx,Ny,Nz)
                inc_diff = Inc_F - Inc_N
            elif abs(field_dir_vec[1]/helpers.norm(field_dir_vec)) ==1: # then allong y
                dH = math.sqrt(Nx**2 + Nz**2) - math.sqrt(MDx**2 + MDz**2)
                dZ = Ny - MDy
                Dec_F, Inc_F, R_F = helpers.cart2dir(field_dir_vec[0],field_dir_vec[1],field_dir_vec[2])
                Dec_N, Inc_N, R_N = helpers.cart2dir(Nx,Ny,Nz)
                inc_diff = Inc_F - Inc_N
            elif abs(field_dir_vec[2]/helpers.norm(field_dir_vec)) ==1: # then allong z
                dH = math.sqrt(Nx**2 + Ny**2) - math.sqrt(MDx**2 + MDy**2)
                dZ = Nz - MDz
                Dec_F, Inc_F, R_F = helpers.cart2dir(field_dir_vec[0],field_dir_vec[1],field_dir_vec[2])
                Dec_N, Inc_N, R_N = helpers.cart2dir(Nx,Ny,Nz)
                inc_diff = Inc_F - Inc_N
            else: # labfiels is not along any of the principal axis, do not calulate dt*
                dH = 9999

            if dH != 9999: # if it is, the labfield is not allong pinciple axix and no calculations are there

                B = dH / math.tan(theta_dt)

                Lim_upp = 2968  # rad *1000 -> 2.968 rad = 170 degrees
                Lim_low = 175   # rad *1000 -> 0.175 rad = 10 degrees

                # check if data is within limits
                # floor rounds to the nearest small integer
                if (math.floor(theta_dt*1000)) < Lim_upp and (math.floor(theta_dt*1000)) > Lim_low:
                    if inc_diff > 0:
                        t_star.append( (-dZ + B) * abs(b_slope) * 100 / abs(Yint))
                    else:
                        t_star.append( (dZ - B) * abs(b_slope) * 100 / abs(Yint))
                if (math.floor(theta_dt*1000)) < Lim_low:
                    t_star.append(0)
                if (math.floor(theta_dt*1000)) > Lim_upp:
                    t_star.append( -dZ * 100 / (abs(Yint) + abs(Xint)))

        if max(t_star)> 0:
            dt_star = max(t_star)
        else:
            dt_star = 0


    sc["check_statistics"]["Tail_check_stat"]["dt_star"] = dt_star
    return sc



def Additivity_check_stat(sc):
    """
    Function to calculate additivity statistics. The number of checks (n_add) and the maximum absolute additivity check difference normalized by the total TRM (obtained from the intersection of the best-fit line and the x-axis on an Arai plot; Leonhardt et al., 2004a).

    input: additivity check info (add_check_steps, add_check_step, add_check_vec), info for the pTRM (x_ptrm_all, x_temp_all), the intercept with the x-asix (Xint)
    output: n_add, d_AC
    """
    # input:    preprocessed/       msrmnts[add_check_steps]
    #                               checks[add_check_step, add_check_vec]
    #                               basics[x_ptrm_all, x_temp_all]
    #           arai_statistics/    intercept_stats[Xint]
    # output:   check_statistics/   Additivity_check[n_add, d_AC]


    # Determine the Additivity check for all checks with Ti < Tmax, that is, all checks that are in the selection, use the ptrm of all steps
    add_check_steps = sc["preprocessed"]["msrmnts"]["add_check_steps"]
    add_check_step = sc["preprocessed"]["checks"]["add_check_step"]
    add_check_vec = sc["preprocessed"]["checks"]["add_check_vec"]

    x_ptrm_all = sc["preprocessed"]["basics"]["x_ptrm_all"]
    x_temp_all = sc["preprocessed"]["basics"]["x_temp_all"]
    Xint = sc["arai_statistics"]["intercept_stats"]["Xint"]

    n_add = []
    d_AC = []

    if len(add_check_steps) != 0:  # then do calculations

        AC = []
        for i in range(len(add_check_vec)):
            for j in range(len(x_ptrm_all)):
                if add_check_step[i][0] == x_temp_all[j]:   # step for the add check and the ptrm should be the same
                    AC.append(helpers.norm(add_check_vec[i]) - x_ptrm_all[j])

        n_add = len(AC)

        AC_abs = []
        for l in range(len(AC)):
            AC_abs.append(abs(AC[l]))

        d_AC = 100* max(AC_abs) / abs(Xint)


    sc["check_statistics"]["Additivity_check"]["n_add"] = n_add
    sc["check_statistics"]["Additivity_check"]["d_AC"] = d_AC
    return sc
