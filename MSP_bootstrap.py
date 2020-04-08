import math
import helpers

"""

This script is used to calculate the Bootstrap parameters. This is independent of the suitcase.

"""


def MSP_boots_calc(site, selection, alpha, NumCycles, Confidence):
    """
    This function is used to get the Bootstrap parameters, and contains two functions: (1)bootstrap and (2) specimen_fail_pass. The first function, bootstrap, is called  for the different options DB, DSC, DSC_corr and returns the bootstrap values: PI_min, PI_max, Boot_int_min, Boot_int_max, Boot_avg.

    The second function is used to check if specimens are outside the bootstrap interval, with the input of Q_str, c_str, site, Boot_min, Boot_max. And returns the specimen name and if it fails or passes.

    input: site selection information, the alpha used for DSC as input by the user, and information for the bootstrap; NumCycles & Confidence
    output: returns the Bootstrap parameters, and the specimen names if they fail or pass for the different options
    """

    # input:    input    [site, selection, alpha, NumCyles, Confidence]
    # output:   B_DB, B_DB_corr, B_DSC, B_DSC_corr, Specimens_DB, Specimens_DB_corr, Specimens_DSC, Specimens_DSC_corr, specimen_fail_pass

    # Do the Bootstrap for both Q_DB & Q_DSC with new function & for the corected values
    def boostrap(Q_str, c_str, site, selection, alpha, NumCycles, Confidence):
        # split in measurements m0 m1 m2 m3 m4 with multiple specimens per list
        m0 = list(filter(lambda m: m['type'] == 0, selection))
        m1 = list(filter(lambda m: m['type'] == 1, selection))
        m2 = list(filter(lambda m: m['type'] == 2, selection))
        m3 = list(filter(lambda m: m['type'] == 3, selection))
        m4 = list(filter(lambda m: m['type'] == 4, selection))

        m1_all = list(filter(lambda m: m['type'] == 1, site))

        # get the steps for the labfield array, this is done by looking at all the data from one site and find the min and maximum used labfields.

        fields = []
        num_specimens = len(m1_all)                   # all the data and not only the selection
        for j in range(num_specimens):
            fields.append(m1_all[j]["lab_field"])   # append all used labfields

        # find min and max labfield used and determine the step
        minField = min(fields)
        maxField = max(fields)
        numsteps = 11                                   # Moster et al., 2015 shows that 11 lab steps give the best results
        step = (minField + maxField) / (numsteps - 1.)    # This is (Hmin+Hmax)/10

        # append the step to a list of labfields -> Hlist
        Hlist = []
        for i in range(numsteps):
            Hlist.append(i * step)

        # set minimum standard deviation for Hlab
        stdevH_min = 10

        N2 = []
        stdevHl = []
        aa = []
        bb = []
        intercept = []

        H0 = []
        H1 = []
        H2 = []
        H3 = []
        H4 = []
        H5 = []
        H6 = []
        H7 = []
        H8 = []
        H9 = []
        H10 = []

        m = 0
        killCounter = 0
        while m < (NumCycles) and killCounter < (NumCycles * 5):

            Hlab_DB = []
            Hlab_DSC = []
            Q_DB_error = []
            Q_DSC_error = []

            num_specimens = len(m0)

            for j in range(num_specimens):  # get N times a random specimen

                # get the index of a random specimen
                i = int(helpers.rand_num() * num_specimens)  # random number between 0 & N

                # get moment per random specimen
                m_m0 = m0[i]["total_m"]
                m_m1 = m1[i]["total_m"]

                # get corresponding error for that specimen, for Q_DB only m0 & m1
                e_m0 = m0[i]["error"]
                e_m1 = m1[i]["error"]

                # calculate new m0_err and m1_err to calculate new Q_DB_error
                frac_m0 = helpers.rand_num() * (0.02 * e_m0) + 1 - 0.01 * e_m0
                m0_err = frac_m0 * m_m0

                frac_m1 = helpers.rand_num() * (0.02 * e_m1) + 1 - 0.01 * e_m1
                m1_err = frac_m1 * m_m1

                Q_DB_error.append((m1_err - m0_err) / m0_err)
                Hlab_DB.append(m1[i]["lab_field"])

                if Q_str == "DSC":
                    if m2[i]["total_m"] != None:

                        m_m2 = m2[i]["total_m"]
                        m_m3 = m3[i]["total_m"]
                        m_m4 = m4[i]["total_m"]

                        e_m2 = m2[i]["error"]
                        e_m3 = m3[i]["error"]
                        e_m4 = m4[i]["error"]

                        # and check for the corrected version, if so replace the moments
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

                        frac_m0 = helpers.rand_num() * (0.02 * e_m0) + 1 - 0.01 * e_m0
                        m0_err = frac_m0 * m_m0

                        frac_m1 = helpers.rand_num() * (0.02 * e_m1) + 1 - 0.01 * e_m1
                        m1_err = frac_m1 * m_m1
                        frac_m2 = helpers.rand_num() * (0.02 * e_m2) + 1 - 0.01 * e_m2
                        m2_err = frac_m2 * m_m2

                        frac_m3 = helpers.rand_num() * (0.02 * e_m3) + 1 - 0.01 * e_m3
                        m3_err = frac_m3 * m_m3

                        Q_DSC_error.append(2 * ((1 + alpha) * m1_err - m0_err - alpha * m3_err) / (2 * m0_err - m1_err - m2_err))
                        Hlab_DSC.append(m2[i]["lab_field"])

            if (Q_str == "DB"):
                Q_error = Q_DB_error
                Hlab = Hlab_DB
            elif (Q_str == "DSC"):
                Q_error = Q_DSC_error
                Hlab = Hlab_DSC

            N = len(Q_error)

            if N > 1:
                avgH = sum(Hlab) / N

                # calculate standard deviation on Hlab, and determine x and y
                stdevH1 = []
                x = []
                y = []
                for k in range(N):
                    stdevH1.append((Hlab[k] - avgH)**2)
                    x.append(Hlab[k])
                    y.append(Q_error[k])
                stdevH = math.sqrt(sum(stdevH1) / (N - 1))

                # calculate Sx, Sy, Sxx, Syy, Sxy
                Sx = sum(x)
                Sy = sum(y)
                Sxy = helpers.dot_product(x, y)
                Sxx = helpers.dot_product(x, x)

                # calculate linear fit is not all at the same Hlab
                if stdevH > stdevH_min:
                    b = (N * Sxy - Sx * Sy) / (N * Sxx - Sx**2)
                    a = Sy / N - b * Sx / N

                    PI = -1 * a / b

                    N2.append(N)
                    stdevHl.append(stdevH)
                    aa.append(a)
                    bb.append(b)
                    intercept.append(PI)

                    H0.append(a + b * Hlist[0])
                    H1.append(a + b * Hlist[1])
                    H2.append(a + b * Hlist[2])
                    H3.append(a + b * Hlist[3])
                    H4.append(a + b * Hlist[4])
                    H5.append(a + b * Hlist[5])
                    H6.append(a + b * Hlist[6])
                    H7.append(a + b * Hlist[7])
                    H8.append(a + b * Hlist[8])
                    H9.append(a + b * Hlist[9])
                    H10.append(a + b * Hlist[10])

                    # end of the big while loop, add one to m (this should be within the if statement)
                    m += 1
            killCounter += 1

        # sort columns and apply cut-off
        cutOffValue = 0.01 * (100 - Confidence) / 2
        cutOff = int(NumCycles * cutOffValue)

        H0.sort()
        H1.sort()
        H2.sort()
        H3.sort()
        H4.sort()
        H5.sort()
        H6.sort()
        H7.sort()
        H8.sort()
        H9.sort()
        H10.sort()

        Q_Hlist = [H0, H1, H2, H3, H4, H5, H6, H7, H8, H9, H10]

        # determine the average of the bootstrap over the 11 labfields
        # take the average of each of the labfield specified in Q_Hlist

        Boot_int_min = []
        Boot_int_max = []
        Boot_avg = []

        if len(Q_Hlist[0]) != 0:
            h = 0
            for el in Q_Hlist:
                Boot_avg.append([Hlist[h], sum(el) / len(el)])
                h += 1

            F = cutOff              # the minimum value F first
            L = m - cutOff - 1      # the maximum value L last ( -1 because python counts from 0)

            y_min = []
            y_max = []
            for w in range(len(Q_Hlist)):
                y_min.append(Q_Hlist[w][F])
                y_max.append(Q_Hlist[w][L])

            for w in range(len(Hlist)):
                Boot_int_min.append([Hlist[w], y_min[w]])
                Boot_int_max.append([Hlist[w], y_max[w]])

            # determine the x axis intercept for lower bound
            ind_min = 999
            for i in range(len(y_min) - 1):
                if (y_min[i] < 0) & (y_min[i + 1] > 0):
                    ind_min = i

            if ind_min == 999:
                ictLow = None
            else:
                slope_min = (y_min[ind_min + 1] - y_min[ind_min]) / (Hlist[ind_min + 1] - Hlist[ind_min])
                ictLow = -1 * (y_min[ind_min] - Hlist[ind_min] * slope_min) / slope_min

            # determine the x axis intercept for upper bound
            ind_max = 999
            for j in range(len(y_max) - 1):
                if (y_max[j] < 0) & (y_max[j + 1] > 0):
                    ind_max = j

            if ind_max == 999:
                ictHigh = None
            else:
                slope_max = (y_max[ind_max + 1] - y_max[ind_max]) / (Hlist[ind_max + 1] - Hlist[ind_max])
                ictHigh = -1 * (y_max[ind_max] - Hlist[ind_max] * slope_max) / slope_max

            # write corresponding PI min and max values, these are the intercepts of the bootstrap intervals
            PI_min = ictHigh
            PI_max = ictLow
        else:
            PI_min = None
            PI_max = None

        return [PI_min, PI_max, Boot_int_min, Boot_int_max, Boot_avg]

    def specimen_fail_pass(Q_str, c_str, site, Boot_min, Boot_max):

        # split in measurements m0 m1 m2 m3 m4 with multiple specimens per list
        m0 = list(filter(lambda m: m['type'] == 0, site))
        m1 = list(filter(lambda m: m['type'] == 1, site))
        m2 = list(filter(lambda m: m['type'] == 2, site))
        m3 = list(filter(lambda m: m['type'] == 3, site))
        m4 = list(filter(lambda m: m['type'] == 4, site))

        # determine if specimen is above / below bootstrap interval and FAIL or PASS specimen

        num_specimens = len(m0)

        # calculate original Q_DB & Q_DSC (also corrected) for ALL specimens
        Hlab = []
        Q_DB = []
        Q_DSC = []
        B_specimen_pass_fail = []
        for i in range(num_specimens):
            # get labfield and name specimen
            Hlab = m1[i]["lab_field"]
            name = m0[i]["specimen"]  # specimen name is the same for m0-m4

            m_m0 = m0[i]["total_m"]
            m_m1 = m1[i]["total_m"]
            m_m2 = m2[i]["total_m"]
            m_m3 = m3[i]["total_m"]
            m_m4 = m4[i]["total_m"]

            if (m_m2 == None) or (m_m3 == None) or (m_m4 == None):
                Q_DB.append([name, Hlab, (m_m1 - m_m0) / m_m0])
                Q_DSC.append([name, None, None])
            else:
                # MSP-DSC calculate  Q_DB & Q_DSC (check for corr)

                Q_DB.append([name, Hlab, (m_m1 - m_m0) / m_m0])  # no Corr

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

                Q_DSC.append([name, Hlab, 2 * ((1 + alpha) * m_m1 - m_m0 - alpha * m_m3) / (2 * m_m0 - m_m1 - m_m2)])

        # loop over all specimens and check closest labfields from the Boot_min and Boot_max
        for i in range(num_specimens):
            Hsam = Q_DB[i][1]
            if (Q_str == "DB"):         # get the y for specimen for Q_DB
                ysam = Q_DB[i][2]
            elif (Q_str == "DSC"):      # get the y for specimen for Q_DSC
                ysam = Q_DSC[i][2]

            name = m0[i]["specimen"]  # specimen name is the same for m0-m4

            if ysam == None:
                B_specimen_pass_fail.append([name, "None"])
            else:
                ind_min = 999
                for j in range(len(Boot_min) - 1):
                    if (Boot_min[j][0] <= Hsam) and (Boot_min[j + 1][0] > Hsam):
                        ind = j

                a_min = (Boot_min[ind + 1][1] - Boot_min[ind][1]) / (Boot_min[ind + 1][0] - Boot_min[ind][0])
                ycalc_min = (Hsam - Boot_min[ind][0]) * a_min + Boot_min[ind][1]

                a_max = (Boot_max[ind + 1][1] - Boot_max[ind][1]) / (Boot_max[ind + 1][0] - Boot_max[ind][0])
                ycalc_max = (Hsam - Boot_max[ind][0]) * a_max + Boot_max[ind][1]

                if (ysam > ycalc_max) or (ysam < ycalc_min):
                    B_specimen_pass_fail.append([name, "fail"])
                else:
                    B_specimen_pass_fail.append([name, "pass"])

        return [B_specimen_pass_fail]

    # run the diffent options for the Bootstrap
    B_DB = boostrap("DB", "", site, selection, alpha, NumCycles, Confidence)
    B_DB_corr = []  # Does not exist

    B_DSC = boostrap("DSC", "", site, selection, alpha, NumCycles, Confidence)
    B_DSC_corr = boostrap("DSC", "_corr", site, selection, alpha, NumCycles, Confidence)

    # find if specimens are outside bootstrap interval with routine above, give the correct upper and lower intervals from the bootstrap
    Specimens_DB = specimen_fail_pass("DB", "", site, B_DB[2], B_DB[3])
    Specimens_DB_corr = []
    Specimens_DSC = specimen_fail_pass("DSC", "", site, B_DSC[2], B_DSC[3])
    Specimens_DSC_corr = specimen_fail_pass("DSC", "_corr", site, B_DSC_corr[2], B_DSC_corr[3])

    return [B_DB, B_DB_corr, B_DSC, B_DSC_corr, Specimens_DB, Specimens_DB_corr, Specimens_DSC, Specimens_DSC_corr]
