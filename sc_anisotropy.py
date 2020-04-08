import math
import numpy
import helpers

"""
This script contains the anisotropy calculations.
"""

def s_tensor_calc(sc):
    """
    Function to calculate the s-tensor for the anisotropy correction, following the calculations in the standard paleointensity definitions by Paterson et al., 2014.

    input: anisotropy measurements
    output: the s-tensor
    """

    # input:    preprocessed/               aniso_trm            [x+, x-, y+, y-, z+, z-, check, orderd]
    # output:   anisotropy_statistics/      s_tensor   [s]

    Xp_list = sc["preprocessed"]["aniso_trm"]["x+"]
    Xm_list = sc["preprocessed"]["aniso_trm"]["x-"]
    Yp_list = sc["preprocessed"]["aniso_trm"]["y+"]
    Ym_list = sc["preprocessed"]["aniso_trm"]["y-"]
    Zp_list = sc["preprocessed"]["aniso_trm"]["z+"]
    Zm_list = sc["preprocessed"]["aniso_trm"]["z-"]
    check_list = sc["preprocessed"]["aniso_trm"]["check"]

    orderd = sc["preprocessed"]["aniso_trm"]["orderd"]


    if Xp_list != None: # then to the calculations for anisotropy tensor
        master_list = [Xp_list, Yp_list, Zp_list, Xm_list, Ym_list, Zm_list]

        # get the order of the aniso measurements and make TRM vec of all measurements
        meas_ord = []
        for meas in orderd:
            if meas["type"] != 87: # do not append the check, type 80 is not in orderd
                meas_ord.append(int(meas["type"])-80)

        # small helper funtion to add coordinates to TRM vec
        def construct_TRM(input_list):
            TRM.append(input_list[0]["x"])
            TRM.append(input_list[0]["y"])
            TRM.append(input_list[0]["z"])
            return TRM

        # add all TRM x, y, z coordinates of the 6 measurements in the correct order
        TRM = []
        for i in range(len(master_list)):
            TRM= construct_TRM(master_list[int(meas_ord[i])-1])

        # make P matrix with the 6 direction positions of measuring the TRM
        xp_dir = [1, 0, 0]      # x+
        yp_dir = [0, 1, 0]      # y+
        zp_dir = [0, 0, 1]      # z+
        xm_dir = [-1, 0, 0]     # x-
        ym_dir = [0, -1, 0]     # y-
        zm_dir = [0, 0, -1]     # z-

        # add all direction to a master list
        master_dir = [xp_dir, yp_dir, zp_dir, xm_dir, ym_dir, zm_dir]

        # make the P matrix with correct directions according to measurement order
        P = []
        for i in range(len(meas_ord)):
            P.append(master_dir[int(meas_ord[i])-1])

        # generate the design matrix A
        A = []
        for i in range(len(P)):
            Ai1 = [P[i][0], 0,                0,    P[i][1],        0,      P[i][2] ]
            Ai2 = [0,       P[i][1],          0,    P[i][0],    P[i][2],    0       ]
            Ai3 = [0,       0,          P[i][2],    0,          P[i][1],    P[i][0] ]
            A.append(Ai1)
            A.append(Ai2)
            A.append(Ai3)

        At = helpers.transpose_list(A)

        AtA = []
        AtA_j = []
        for i in range(len(At)):
            for j in range(len(At)):
                AtA_j.append(helpers.dot_product(At[i],At[j]))
            AtA.append(AtA_j)
            AtA_j = []

        # make numpy arrays and do numpy linalg.inv
        # make numpy array (1)
        AtA = numpy.array(AtA)

        # do numpy calculation (2)
        AtA_inv= numpy.linalg.inv(AtA)


        # back to python lists (3)
        AtA_inv = AtA_inv.tolist()

        # make x, that is the inversion of AtA times At
        x = []
        x_j = []
        for i in range(len(AtA_inv)):
            for j in range(len(A)):
                x_j.append(helpers.dot_product(AtA_inv[i],A[j]))
            x.append(x_j)
            x_j = []

        # next step is to muliply x with TRM tensor to get s
        # s (6x1) = x (6 x 18) * TRM (18 x 1)
        s = []
        for i in range(len(x)):
            s.append(helpers.dot_product(x[i], TRM))

        sc["anisotropy_statistics"]["Aniso_tensor"]["s_tensor"] = s
    return(sc)


def anisotropy_calc(sc):
    """
    Function that calculates the paleointensity correction factor when anisotropy data is available. It uses the previously calculated s-tensor. First the direction of the NRM for the selected part of the Arai plot is determined to give Mhat_ChRM. Which is used to get the anisotropy correction c. The anisotropy correction in multiplied with Banc to get Banc_aniso_corr.

    input: s-tensor, the direction of the applied labfield, direction of the data selection, the original Banc
    output: the anisotropy correction anis_c, the anisotropy corrected paleointensity estimate Banc_aniso_corr
    """

    # input:    preprocessed/               field_basics            [field_dir_vec]
    #                                       s_tensor                [s1_list, s2_list, s3_list, s4_list, s5_list, s6_list, scheck_list]
    #           arai_statistics/            PI_Banc_est             [B_anc]
    #           directional_statistics/     mean_dir_stat           [Mdec_free, Minc_free]
    # output:   anisotropy_statistics/      Anisotropy_Correction   [anis_c, Banc_aniso_corr]

    Mdec_free = sc["directional_statistics"]["mean_dir_stat"]["Mdec_free"]
    Minc_free = sc["directional_statistics"]["mean_dir_stat"]["Minc_free"]

    field_dir_vec = sc["preprocessed"]["field_basics"]["field_dir_vec"]
    B_anc = sc["arai_statistics"]["PI_Banc_est"]["B_anc"]

    Xp_list = sc["preprocessed"]["aniso_trm"]["x+"]
    s_tensor = sc["anisotropy_statistics"]["Aniso_tensor"]["s_tensor"]

    if Xp_list != None: # then there is anisotropy data, do calculations, if no anisotropy data present do nothing

        anis_c = []

        # the direction of the NRM, should be determined from the selected Arai plot
        # this is the Mdec_free & Minc_free, the unit vecotr of this gives -> Mhat_ChRM

        # get mean unit vector
        Mhat_ChRM = helpers.dir2cart(Mdec_free, Minc_free, 1)
        Blab_orient = field_dir_vec # it was: helpers.dir2cart(Mdec_free, Minc_free, 1), that is wrong!! (13 feb 2020)

        s1 = s_tensor[0]
        s2 = s_tensor[1]
        s3 = s_tensor[2]
        s4 = s_tensor[3]
        s5 = s_tensor[4]
        s6 = s_tensor[5]

        A1 = [s1, s4, s6]
        A2 = [s4, s2, s5]
        A3 = [s6, s5, s3]

        A = [A1, A2, A3]

        # make A and Mhat_ChRM into a numpy arrays (1)
        A = numpy.array(A)
        Mhat_ChRM = numpy.array(Mhat_ChRM)

        # do numpy calculation (2)
        Hanc = numpy.linalg.solve(A, Mhat_ChRM)

        # back to python lists (3)
        A = A.tolist()
        Mhat_ChRM = Mhat_ChRM.tolist()
        Hanc = Hanc.tolist()

        # unit vector in the direction of the ancient field
        Hanc_hat = helpers.list_div_num(Hanc, helpers.norm(Hanc))

        Manc = []
        Mlab = []
        for i in range(len(A)):
            Manc.append(helpers.dot_product(A[i],Hanc_hat))
            Mlab.append(helpers.dot_product(A[i],field_dir_vec))

        aniso_c = helpers.norm(Mlab)/helpers.norm(Manc)

        Banc_aniso_corr = aniso_c*B_anc

        sc["anisotropy_statistics"]["Anisotropy_Correction"]["aniso_c"] = aniso_c
        sc["anisotropy_statistics"]["Anisotropy_Correction"]["Banc_aniso_corr"] = Banc_aniso_corr

    return(sc)

def aniso_alteration_calc(sc):
    """
    Function that calculates the anisotropy alteration check, dTRM_anis.

    input: anisotropy measurements and the order they are measured
    output: the anisotropy alteration dTRM_anis
    """

    # input:    preprocessed/aniso_trm                  [x+, x-, y+, y-, z+, z-, check, orderd]
    #           input/                                  [specimen]
    #
    # output:   anisotropy_statistics/aniso_alteration  [dTRM_anis]

    Xp_list = sc["preprocessed"]["aniso_trm"]["x+"]
    Yp_list = sc["preprocessed"]["aniso_trm"]["y+"]
    Zp_list = sc["preprocessed"]["aniso_trm"]["z+"]
    Xm_list = sc["preprocessed"]["aniso_trm"]["x-"]
    Ym_list = sc["preprocessed"]["aniso_trm"]["y-"]
    Zm_list = sc["preprocessed"]["aniso_trm"]["z-"]
    check_list = sc["preprocessed"]["aniso_trm"]["check"]

    master_list = [Xp_list, Yp_list, Zp_list, Xm_list, Ym_list, Zm_list]

    orderd = sc["preprocessed"]["aniso_trm"]["orderd"]


    specimen = sc["input"]["specimen"]

    dTRM_anis = []

    # check for the direction of the field for the check step and see which step is the same

    if check_list != None:
        if len(check_list) != 0: # only when a check is measured, calculate the pTRM_anis

            # find the used lab directions in the order they are measured
            meas_lab_dir = []
            meas_orderd = []
            for meas in orderd:
                meas_lab_dir.append([int(meas["lab_field_dec"]),int(meas["lab_field_inc"])])
                meas_orderd.append(int(meas["type"]))
            # the check step should be the last measured lab direction, check this!
            print("meas_orderd",meas_orderd)


            if meas_orderd[len(meas_orderd)-1] == 87:

                labdir_check = meas_lab_dir[len(meas_lab_dir)-1]

                # now find the corresponding step, this is when the labfield is the same as used for the check step
                for i in range(len(meas_lab_dir)-1): # -1 because you do not want to take the step into account
                    if meas_lab_dir[i][0] == labdir_check[0] and meas_lab_dir[i][1] == labdir_check[1]:
                        P1 = orderd[i]


                # do the calculation with the total magnetization of P1 and the check step
                dTRM_anis = (abs(P1["total_m"] - check_list[0]["total_m"]) / (P1["total_m"])) *100

        sc["anisotropy_statistics"]["aniso_alteration"]["dTRM_anis"] = dTRM_anis
    return(sc)
