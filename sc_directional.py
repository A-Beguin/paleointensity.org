import math
import numpy
import helpers

"""
This script contains the functions for the directional statistics, according to the SPD (Standard Paleointensity Definitions).

Paterson, G. A., L. Tauxe, A. J. Biggin, R. Shaar, and L. C. Jonestrask (2014), On improving the selection of Thellier-type paleointensity data, Geochem. Geophys. Geosyst., doi: 10.1002/2013GC005135
"""


def mean_dir_stat(sc):
    """
    Function to calculate the direction statistics from the NRM remaining measuremtens.
    This function includes another function "direction_stat_dir_type" this is calles to obtain the Dec, Inc, MAD for the free floating and anchored fit.

    input: NRM remaining
    output: Mdec_free, Minc_free, MAD_free, Mdec_anc, Minc_anc, MAD_anc, x, y, z, CMvec
    """

    # input:    preprocessed/           basics[NRM_rem]
    # output:   directional_statistics/ mean_dir_stat[Mdec_free, Minc_free, MAD_free, Mdec_anc, Minc_anc, MAD_anc, x, y, z, CMvec]

    NRM_rem = sc["preprocessed"]["basics"]["NRM_rem"]

    # print(NRM_rem)
    # get the NRM x, y and z values
    x = []
    y = []
    z = []
    for step in NRM_rem:
        x.append(step["x"])
        y.append(step["y"])
        z.append(step["z"])

    # calculations to get dec inc and MAD, type_dir specifies the free or anchored fit
    def direction_stat_dir_type(x, y, z, type_dir):

        # center of mass
        mean_x = sum(x) / len(x)
        mean_y = sum(y) / len(y)
        mean_z = sum(z) / len(z)

        # transform the NRM vector
        if type_dir == 1:  # free

            x_prime = helpers.list_min_num(x, mean_x)
            y_prime = helpers.list_min_num(y, mean_y)
            z_prime = helpers.list_min_num(z, mean_z)

        if type_dir == 2:  # forced / anchored
            x_prime = x
            y_prime = y
            z_prime = z

        orient_tensor = [[helpers.dot_product(x_prime, x_prime), helpers.dot_product(x_prime, y_prime), helpers.dot_product(x_prime, z_prime)],
                         [helpers.dot_product(x_prime, y_prime), helpers.dot_product(y_prime, y_prime), helpers.dot_product(y_prime, z_prime)],
                         [helpers.dot_product(x_prime, z_prime), helpers.dot_product(z_prime, y_prime), helpers.dot_product(z_prime, z_prime)]]

        # this numpy statement stays here :-)
        orient_tensor = numpy.array(orient_tensor)

        # get the eigenvalues (tau), and eigenvectors(V)
        tau, V = numpy.linalg.eig(orient_tensor)

        # tau & V to list, instead of numpy
        tau = tau.tolist()
        V = V.tolist()

        # print("help",tau, V)

        # rescale tau to sum-to-one
        tau = helpers.list_div_num(tau, sum(tau))

        # find index max tau (eigenvaule 1)
        ind_tau_max = -1
        tau_max = -1
        for idx in range(len(tau)):
            if tau[idx] > tau_max:
                tau_max = tau[idx]
                ind_tau_max = idx

        # transpose V
        TV = helpers.transpose_list(V)

        # find eigenvector and eigenvalue 1
        v1 = TV[ind_tau_max]  # [0] is necessary to make it a vector instead of matrix
        e1 = tau[ind_tau_max]

        # the other two eigenvaules and vectors
        sum_e23 = sum(tau) - tau_max

        # define the reference vector
        r1 = (x_prime[0] - x_prime[len(x_prime) - 1])
        r2 = (y_prime[0] - y_prime[len(x_prime) - 1])
        r3 = (z_prime[0] - z_prime[len(x_prime) - 1])

        R = [r1, r2, r3]

        dot = helpers.dot_product(v1, R)

        if dot < -1:
            dot = -1
        elif dot > 1:
            dot = 1

        if math.acos(dot) > (math.pi / 2.):
            PD = helpers.list_mult_num(v1, -1)
        else:
            PD = v1

        [Mdec, Minc, R] = helpers.cart2dir(PD[0], PD[1], PD[2])

        CMvec = [mean_x, mean_y, mean_z]

        MAD = math.atan(math.sqrt(sum_e23 / e1)) * 180 / math.pi
        return [Mdec, Minc, MAD, CMvec]

    [Mdec_free, Minc_free, MAD_free, CMvec] = direction_stat_dir_type(x, y, z, 1)
    [Mdec_anc, Minc_anc, MAD_anc, CMvec_anc] = direction_stat_dir_type(x, y, z, 2)

    sc["directional_statistics"]["mean_dir_stat"]["x"] = x
    sc["directional_statistics"]["mean_dir_stat"]["y"] = y
    sc["directional_statistics"]["mean_dir_stat"]["z"] = z
    sc["directional_statistics"]["mean_dir_stat"]["Mdec_free"] = Mdec_free
    sc["directional_statistics"]["mean_dir_stat"]["Minc_free"] = Minc_free
    sc["directional_statistics"]["mean_dir_stat"]["MAD_free"] = MAD_free
    sc["directional_statistics"]["mean_dir_stat"]["Mdec_anc"] = Mdec_anc
    sc["directional_statistics"]["mean_dir_stat"]["Minc_anc"] = Minc_anc
    sc["directional_statistics"]["mean_dir_stat"]["MAD_anc"] = MAD_anc
    sc["directional_statistics"]["mean_dir_stat"]["CMvec"] = CMvec
    return sc


def NRM_dir_stat(sc):
    """
    Function to calculate the NRM direction statistics, DANG, free floating direction and NRMdev.

    input: Mdec_free,Minc_free, x, y, z, Yint
    output: DANG, NRMdev, dir_vec_free
    """

    # input:    directional_statistics/ mean_dir_stat[Mdec_free,Minc_free, x, y, z]
    #           arai_statistics/        intercept_stats[Yint]
    # output:   directional_statistics/ NRM_dir_stat[DANG, NRMdev, dir_vec_free]
    x = sc["directional_statistics"]["mean_dir_stat"]["x"]
    y = sc["directional_statistics"]["mean_dir_stat"]["y"]
    z = sc["directional_statistics"]["mean_dir_stat"]["z"]
    Mdec_free = sc["directional_statistics"]["mean_dir_stat"]["Mdec_free"]
    Minc_free = sc["directional_statistics"]["mean_dir_stat"]["Minc_free"]
    Yint = sc["arai_statistics"]["intercept_stats"]["Yint"]

    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    mean_z = sum(z) / len(z)

    center_of_mass = [mean_x, mean_y, mean_z]

    # calculate DANG & NRMdev
    dir_vec_free = helpers.dir2cart(Mdec_free, Minc_free, 1)

    DANG = helpers.get_angle_diff(dir_vec_free, center_of_mass)

    NRMdev = ((math.sin(math.radians(DANG)) * helpers.norm(center_of_mass)) / abs(Yint)) * 100

    sc["directional_statistics"]["NRM_dir_stat"]["dir_vec_free"] = dir_vec_free
    sc["directional_statistics"]["NRM_dir_stat"]["DANG"] = DANG
    sc["directional_statistics"]["NRM_dir_stat"]["NRMdev"] = NRMdev
    return sc


def alpha_stat(sc):
    """
    Function to calculate the alpha statistic, the angular difference between the anchored and free-floating best-fit directions

    input: free floating direction vector (dir_vec_free), and the anchored direction (Mdec_anc, Minc_anc)
    output: alpha
    """

    # input:    directional_statistics/ mean_dir_stat[Mdec_anc, Minc_anc]
    #                                   NRM_dir_stat[dir_vec_free]
    # output:   directional_statistics/ alpha_stat[alpha]
    Mdec_anc = sc["directional_statistics"]["mean_dir_stat"]["Mdec_anc"]
    Minc_anc = sc["directional_statistics"]["mean_dir_stat"]["Minc_anc"]
    dir_vec_free = sc["directional_statistics"]["NRM_dir_stat"]["dir_vec_free"]

    dir_vec_anc = []
    dir_vec_anc = helpers.dir2cart(Mdec_anc, Minc_anc, 1)

    alpha = helpers.get_angle_diff(dir_vec_free, dir_vec_anc)

    sc["directional_statistics"]["alpha_stat"]["alpha"] = alpha
    return sc


def theta_stat(sc):
    """
    Function to calculate the theta statistic, the angle between the applied field direction and the ChRM direction of the NRM

    input: free floating direction of best-fit selection (dir_vec_free), and the applied field direction vector (field_dir_vec)
    output: the angle between the two (theta)
    """

    # input:    directional_statistics/ NRM_dir_stat[dir_vec_free]
    #           preprocessed/           field_basics[field_dir_vec]
    # output:   directional_statistics/ theta_stat[theta]
    dir_vec_free = sc["directional_statistics"]["NRM_dir_stat"]["dir_vec_free"]
    field_dir_vec = sc["preprocessed"]["field_basics"]["field_dir_vec"]

    theta = helpers.get_angle_diff(dir_vec_free, field_dir_vec)

    sc["directional_statistics"]["theta_stat"]["theta"] = theta
    return sc


def gamma_stat(sc):
    """
    Function to calculate the gamma statistic, The angle between the pTRM acquired at the last step used for the best-fit segment and the applied field direction (BLab)

    input: pTRM vector (ptrm), and the applied field direction vector (field_dir_vec)
    output: paleointensity estimate (B_anc)
    """
    # input:    preprocessed/           basics[ptrm]
    #                                   field_basics[field_dir_vec]
    # output:   directional_statistics/ gamma_stat[gamma]
    ptrm = sc["preprocessed"]["basics"]["ptrm"]
    field_dir_vec = sc["preprocessed"]["field_basics"]["field_dir_vec"]

    ptrm_last_step = []
    ptrm_last_step = [ptrm[len(ptrm) - 1]["x"], ptrm[len(ptrm) - 1]["y"], ptrm[len(ptrm) - 1]["z"]]

    gamma = helpers.get_angle_diff(ptrm_last_step, field_dir_vec)

    sc["directional_statistics"]["gamma_stat"]["gamma"] = gamma
    return sc
