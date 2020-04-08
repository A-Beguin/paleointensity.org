import math
import helpers

"""
In this script the input data will be split on the type of measurements. The NRM and pTRM are calculated, as are other basic parameters and field basics that are used in multiple functions in different scripts.
"""

def splitup_msrmnts(sc):
    """
    Function to split-up the input in measurements (msrmts), this is done by looking at the type of measurement and the selection or entire specimen.

    input: input tables from the website, or the demo data. Input from specimen and the selection of the specimen.
    output: the measurements split on type, e.g. zerofield, infield etc.
    """

    # input:  input             [specimen, specimen_selection]
    # output: preprocessed/     msrmnts[zerofield_steps, zerofield_steps_all, infield_steps, infield_antiparallel_steps, ptrm_check_steps,
    #                                   tail_check_steps, ARM_acq_steps, ARM_acq_steps_all, ARM_dem_steps]

    specimen = sc["input"]["specimen"]
    specimen_selection = sc["input"]["specimen_selection"]

    # Split up the different tyes of measurements
    sc["preprocessed"]["msrmnts"]["zerofield_steps"] = list(filter(lambda m: m['type'] == 0, specimen_selection))
    sc["preprocessed"]["msrmnts"]["zerofield_steps_all"] = list(filter(lambda m: m['type'] == 0, specimen))

    sc["preprocessed"]["msrmnts"]["infield_steps"] = list(filter(lambda m: m['type'] == 1, specimen_selection))
    sc["preprocessed"]["msrmnts"]["infield_steps_all"] = list(filter(lambda m: m['type'] == 1, specimen))

    sc["preprocessed"]["msrmnts"]["ptrm_check_steps"] = list(filter(lambda m: m['type'] == 2, specimen_selection))
    sc["preprocessed"]["msrmnts"]["tail_check_steps"] = list(filter(lambda m: m['type'] == 3, specimen_selection))
    sc["preprocessed"]["msrmnts"]["add_check_steps"] = list(filter(lambda m: m['type'] == 4, specimen_selection))

    # Thellier-Thellier
    sc["preprocessed"]["msrmnts"]["infield_antiparallel_steps"] = list(filter(lambda m: m['type'] == 5, specimen_selection))
    sc["preprocessed"]["msrmnts"]["infield_antiparallel_steps_all"] = list(filter(lambda m: m['type'] == 5, specimen))

    # pseudo Thellier
    sc["preprocessed"]["msrmnts"]["ARM_acq_steps"] = list(filter(lambda m: m['type'] == 6, specimen_selection))
    sc["preprocessed"]["msrmnts"]["ARM_acq_steps_all"] = list(filter(lambda m: m['type'] == 6, specimen))
    sc["preprocessed"]["msrmnts"]["ARM_dem_steps"] = list(filter(lambda m: m['type'] == 7, specimen_selection))
    return sc

def prep_anisotropy_tensor(sc):
    """
    Function to prepare the data for the anisotropy correction. Look though the specimen to find the anisotropy types (81 - 87), and place the measurements in the correct arrays of x+, x-, y+, etc.

    input: Specimen data
    output: the anisotropy trm measurements and the order they appear in the input
    """
    # input:  input             [specimen]
    # output: preprocessed/     aniso_trm[orderd, x+, x-, y+, y-, z+, z-, check]

    specimen = sc["input"]["specimen"]

    # get the wanted steps for the anisotropy correction
    aniso_types = [81,82,83,84,85,86,87]
    wanted_types =  aniso_types + [1,5]

    aniso_step = None
    for msrmt in specimen:
        if int(msrmt['type']) in aniso_types:
            aniso_step = msrmt['step']


    # write the ordered input to the suitcase
    orderd = list(filter(lambda m: m['step'] == aniso_step and int(m['type']) in wanted_types, specimen))

    # the anisotropy measurments
    Xp = list(filter(lambda s: int(s['type']) == 81, orderd))
    Yp = list(filter(lambda s: int(s['type']) == 82, orderd))
    Zp = list(filter(lambda s: int(s['type']) == 83, orderd))
    Xm = list(filter(lambda s: int(s['type']) == 84, orderd))
    Ym = list(filter(lambda s: int(s['type']) == 85, orderd))
    Zm = list(filter(lambda s: int(s['type']) == 86, orderd))

    # form the above 6 measuremnts two could miss, these are the field step from the experiment
    p = list(filter(lambda s: int(s['type']) == 1, orderd))
    m = list(filter(lambda s: int(s['type']) == 5, orderd))

    if (len(Xp) + len(Yp) + len(Zp)) != 0: # then you have anisotropy data, write to suitcase

        # check in which direction field was applied and change the 1, 5 for that direction to 8x and 8x
        # check lengths of the anisotropy measurements and put extra information if necessary
        if len(Xm) == 0:
            Xp = p
            Xm = m
            Xp[0]["type"] = 81
            Xm[0]["type"] = 84
            for msrmt in orderd:
                if msrmt["type"] == 1:
                    msrmt["type"] = 81
                elif msrmt["type"] == 5:
                    msrmt["type"] = 84
        elif len(Ym) == 0:
            Yp = p
            Ym = m
            Yp[0]["type"] = 82
            Ym[0]["type"] = 85
            for msrmt in orderd:
                if msrmt["type"] == 1:
                    msrmt["type"] = 82
                elif msrmt["type"] == 5:
                    msrmt["type"] = 85
        elif len(Zm) == 0:
            Zp = p
            Zm = m
            Zp[0]["type"] = 83
            Zm[0]["type"] = 86
            for msrmt in orderd:
                if msrmt["type"] == 1:
                    msrmt["type"] = 83
                elif msrmt["type"] == 5:
                    msrmt["type"] = 86

        orderd = list(filter(lambda s: int(s['type']) != 1 and int(s['type']) != 5, orderd))

        sc["preprocessed"]["aniso_trm"]["x+"] = Xp
        sc["preprocessed"]["aniso_trm"]["x-"] = Xm
        sc["preprocessed"]["aniso_trm"]["y+"] = Yp
        sc["preprocessed"]["aniso_trm"]["y-"] = Ym
        sc["preprocessed"]["aniso_trm"]["z+"] = Zp
        sc["preprocessed"]["aniso_trm"]["z-"] = Zm
        sc["preprocessed"]["aniso_trm"]["check"] = list(filter(lambda s: int(s['type']) == 87, orderd))

        sc["preprocessed"]["aniso_trm"]["orderd"] = orderd
    return(sc)


def calc_nrm_ptrm(sc):
    """
    Function to calculate/split the basic nrm and ptrm vectors and measurements for the different methods and protocols.

    input: preprocessed measurement data, for all protocols and methods from the suitcase
    output: basic vectors and list, like the NRM remaining and the ptrm vectors
    """

    # input:    preprocessed/       msrmnts[zerofield_steps, zerofield_steps_all, infield_steps, infield_steps_all, infield_antiparallel_steps, infield_antiparallel_steps_all, ARM_acq_steps, ARM_acq_steps_all]
    # output:   preprocessed/       basics[nrm0, ptrm, NRM_rem, NRM_vec_all, ptrm_vec_all, x_ptrm_all, x_temp_all, y_nrm_all, y_temp_all, pstep_all]

    zerofield_steps = sc["preprocessed"]["msrmnts"]["zerofield_steps"]
    zerofield_steps_all = sc["preprocessed"]["msrmnts"]["zerofield_steps_all"]
    infield_steps = sc["preprocessed"]["msrmnts"]["infield_steps"]

    infield_steps_all = sc["preprocessed"]["msrmnts"]["infield_steps_all"]
    infield_antiparallel_steps = sc["preprocessed"]["msrmnts"]["infield_antiparallel_steps"]    # Thellier - Thellier
    infield_antiparallel_steps_all = sc["preprocessed"]["msrmnts"]["infield_antiparallel_steps_all"]    # Thellier - Thellier

    ARM_acq_steps = sc["preprocessed"]["msrmnts"]["ARM_acq_steps"]      # pseudo Thellier NAA format
    ARM_acq_steps_all = sc["preprocessed"]["msrmnts"]["ARM_acq_steps_all"]      # pseudo Thellier NAA format

    # get NRM0, the first zerofield step, this is true for both modified Thellier and classical TT and pseudo Thellier (GF & NAA)
    nrm0 = zerofield_steps_all[0]["total_m"]


    NRMrem_list = []
    ptrm_list = []
    NRM_vec_all = []
    ptrm_vec_all = []
    pstep_all = []
    x_ptrm_all = []
    x_temp_all = []
    y_nrm_all = []
    y_temp_all = []


    # "classical thellier Thelier-Thellier (TT)"
    if len(infield_antiparallel_steps) != 0:
        # first the NRM, which is always a zerofield step, also for the Thellier - Thellier protocol,
        # but not always in the selection so check the length to see if it is in the selection

        if len(zerofield_steps) != 0: # then in selection
            NRMrem_list.append([zerofield_steps[0]["step"], zerofield_steps[0]["x"], zerofield_steps[0]["y"], zerofield_steps[0]["z"]])
            ptrm_list.append([zerofield_steps[0]["step"], 0, 0, 0]) # this is new 06-02

        if len(zerofield_steps_all) != 0: # this is probably always the case # this is new 06-02
            NRM_vec_all.append([zerofield_steps_all[0]["x"], zerofield_steps_all[0]["y"], zerofield_steps_all[0]["z"]])
            ptrm_vec_all.append([0,0,0])
            x_ptrm_all.append(0)
            y_nrm_all.append(helpers.norm([zerofield_steps_all[0]["x"], zerofield_steps_all[0]["y"], zerofield_steps_all[0]["z"]]))
            x_temp_all.append(zerofield_steps_all[0]["step"])
            y_temp_all.append(zerofield_steps_all[0]["step"])

        # determine the NRM remaining & ptrm gained for the selection
        for i_step in infield_steps:
            for ai_step in infield_antiparallel_steps:
                if i_step["step"] == ai_step["step"]:
                    NRMrem_list.append([i_step["step"], (i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.])
                    ptrm_list.append([i_step["step"], (i_step["x"] - ai_step["x"]) / 2., (i_step["y"] - ai_step["y"]) / 2., (i_step["z"] - ai_step["z"]) / 2.])


        # determine the NRM remaining & ptrm gained for ALL steps
        for i_step in infield_steps_all:
            for ai_step in infield_antiparallel_steps_all:
                if i_step["step"] == ai_step["step"]:

                    tempor_ptrm = [(i_step["x"] - ai_step["x"]) / 2., (i_step["y"] - ai_step["y"]) / 2., (i_step["z"] - ai_step["z"]) / 2.]
                    tempor_nrm = [(i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.]

                    NRM_vec_all.append(tempor_nrm)
                    ptrm_vec_all.append(tempor_ptrm)
                    x_ptrm_all.append(helpers.norm(tempor_ptrm))
                    x_temp_all.append(i_step["step"])
                    y_nrm_all.append(helpers.norm(tempor_nrm))
                    y_temp_all.append(i_step["step"])


    # "modified thellier"
    if len(infield_antiparallel_steps) == 0:

        # pseudo Thellier GF or Thermal Thellier or Microwave Thellier
        if len(ARM_acq_steps) == 0:
            # if there are more zerofiels steps that infield steps, add steps to the ptrm_list with zeros, assumption is that extra zerofield steps are always before first ptrm step
            # for i in range(len(zerofield_steps) - len(infield_steps)):

            if (len(zerofield_steps) - len(infield_steps)) > 0:
                ptrm_list.append([zerofield_steps_all[0]["step"], 0, 0, 0])  # this takes care off the first y-axis point

            # make ptrm_list
            for i_step in infield_steps:
                for z_step in zerofield_steps:
                    if i_step["step"] == z_step["step"]:
                        ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

        # pseudo-Thellier NAA
        if len(ARM_acq_steps) != 0:
            # get the first ARM step, this is your NRM remaining and the ptrm should be corrected for this
            x0 = ARM_acq_steps_all[0]["x"]
            y0 = ARM_acq_steps_all[0]["y"]
            z0 = ARM_acq_steps_all[0]["z"]
            for a_step in ARM_acq_steps:
                ptrm_list.append([a_step["step"], a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0])

        # make NRMrem_list
        for z_step in zerofield_steps:
            NRMrem_list.append([z_step["step"], z_step["x"], z_step["y"], z_step["z"]])

        for all_step in zerofield_steps_all:
            NRM_vec_all.append([all_step["x"], all_step["y"], all_step["z"]])
            y_nrm_all.append(helpers.norm([all_step["x"], all_step["y"], all_step["z"]]))
            y_temp_all.append(all_step["step"])
            for all_infield in infield_steps_all:       # pTh-NAA is not getting through this if statement
                if all_infield["step"] == all_step["step"]:
                    ptrm_vec_all.append([all_infield["x"] - all_step["x"] ,all_infield["y"] - all_step["y"], all_infield["z"] - all_step["z"]])
                    pstep_all.append(all_infield["step"])
                    x_temp_all.append(all_infield["step"])
                    x_ptrm_all.append(helpers.norm([all_infield["x"] - all_step["x"] ,all_infield["y"] - all_step["y"], all_infield["z"] - all_step["z"]]))

        for a_step in ARM_acq_steps_all:    # pTh -NAA format ptrm vec
            ptrm_vec_all.append([a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0])
            x_ptrm_all.append(helpers.norm([a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0]))
            x_temp_all.append(a_step["step"])

    # rewrite format
    ptrm = []
    for msrmt in ptrm_list:
        ptrm.append({
            "step":msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })

    NRM_rem = []
    for msrmt in NRMrem_list:
        NRM_rem.append({
            "step": msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })

    sc["preprocessed"]["basics"]["nrm0"] = nrm0
    sc["preprocessed"]["basics"]["ptrm"] = ptrm
    sc["preprocessed"]["basics"]["NRM_rem"] = NRM_rem
    sc["preprocessed"]["basics"]["NRM_vec_all"] = NRM_vec_all
    sc["preprocessed"]["basics"]["ptrm_vec_all"] = ptrm_vec_all
    sc["preprocessed"]["basics"]["x_ptrm_all"] = x_ptrm_all
    sc["preprocessed"]["basics"]["x_temp_all"] = x_temp_all
    sc["preprocessed"]["basics"]["y_nrm_all"] = y_nrm_all
    sc["preprocessed"]["basics"]["y_temp_all"] = y_temp_all
    sc["preprocessed"]["basics"]["pstep_all"] = pstep_all
    return sc

def basics(sc):
    """
    Function to calculate more basics from the ptrm and NRM_rem, like the  nrm and ptrm vectors and lists for the different methods and protocols. Vectors are list wiht "x", "y", and "z" values, list are the moments of the vecotr, e.g. x_ptrm_all is a list with all ptrm values for that specimen the "x" stands for the x-axis

    input: the ptrm and NRM_rem vectors with the "step", "x", "y", and "z" values
    output: basic vectors and list, also the average x and y (xBar, yBar)
    """

    # input:    preprocessed/    basics[ptrm, NRM_rem]
    # output:   preprocessed/    basics[y_nrm, xy_temp, NRM_vec_select, x_ptrm, x_temp, ptrm_gained_vec, xBar, yBar]

    ptrm = sc["preprocessed"]["basics"]["ptrm"]
    NRM_rem = sc["preprocessed"]["basics"]["NRM_rem"]


    y_nrm = []
    xy_temp = []
    NRM_vec_select = []
    for nstep in NRM_rem:
        for pstep in ptrm:
            if nstep["step"] == pstep["step"]:
                y_nrm.append(math.sqrt(nstep["x"]**2 + nstep["y"]**2 + nstep["z"]**2))
                xy_temp.append(nstep["step"])
                NRM_vec_select.append([nstep["x"], nstep["y"], nstep["z"]])



    x_ptrm = []
    x_temp = []
    ptrm_gained_vec = []
    for step in ptrm:
        x_ptrm.append(math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2))
        x_temp.append(step["step"])
        ptrm_gained_vec.append([step["x"], step["y"], step["z"]])  # the ptrm gained vector with [x,y,z] per (temp/MW) step

    # calculate average x, and average y -> xBar , yBar
    xBar = sum(x_ptrm) / len(x_ptrm)
    yBar = sum(y_nrm) / len(y_nrm)


    sc["preprocessed"]["basics"]["y_nrm"] = y_nrm
    sc["preprocessed"]["basics"]["xy_temp"] = xy_temp
    sc["preprocessed"]["basics"]["NRM_vec_select"] = NRM_vec_select
    sc["preprocessed"]["basics"]["x_ptrm"] = x_ptrm
    sc["preprocessed"]["basics"]["x_temp"] = x_temp
    sc["preprocessed"]["basics"]["ptrm_gained_vec"] = ptrm_gained_vec
    sc["preprocessed"]["basics"]["xBar"] = xBar
    sc["preprocessed"]["basics"]["yBar"] = yBar
    return sc

def basics_pTh(sc):
    """
    Function to calculate more basics from the ARM_dem_steps and ARM_acq_steps_all for the pseudo-Thellier method [NAA protocol only, NRMdemag, ARM-acquisition, ARMdemag]. This addes on the previous function that calculates the NRM demag basics.

    input: the ARM_dem_steps and ARM_acq_steps_all vectors with the "step", "x", "y", and "z" values
    output: basic vectors and lists for the ARM demag
    """
    # input:    preprocessed/   msrmnts[ARM_dem_steps, ARM_acq_steps_all]
    # output:   preprocessed/   basics_pTh[arm_step, y_arm_d, ARMd_vec,yBar_arm ]

    ARM_dem_steps = sc["preprocessed"]["msrmnts"]["ARM_dem_steps"]      # pseudo Thellier NAA format
    ARM_acq_steps_all = sc["preprocessed"]["msrmnts"]["ARM_acq_steps_all"]

    if len(ARM_dem_steps) != 0: # pseudo-Thellier NAA
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]

        y_arm_d = []
        ARMd_vec = []
        arm_step = []
        for dstep in ARM_dem_steps:
            ARMd_vec.append([dstep["x"]- x0, dstep["y"]- y0, dstep["z"]- z0])
            y_arm_d.append(helpers.norm([dstep["x"] - x0, dstep["y"] -y0 , dstep["z"] - z0]))
            arm_step.append(dstep["step"])


        # calculate average
        yBar_arm = sum(y_arm_d) / len(y_arm_d)

        sc["preprocessed"]["basics_pTh"]["arm_step"] = arm_step
        sc["preprocessed"]["basics_pTh"]["y_arm_d"] = y_arm_d
        sc["preprocessed"]["basics_pTh"]["ARMd_vec"] = ARMd_vec
        sc["preprocessed"]["basics_pTh"]["yBar_arm"] = yBar_arm
    return sc


def field_basics(sc):
    """
    Function to write the field basics to the suitcase, like the lab field strength (Blab), and the direction of the labfield as provided by the user in the input (field_dir_vec_initial). The labfield direction is checked with the direction of the last ptrm_step to check if the field_dir_vec_initial is correct or if it should be anti-parallel to that.

    input: the ptrm_gained_vec from the basics and the infield_steps from the measurements
    output: the field basics: Blab, field_dir_vec_initial, field_dir_vec
    """

    # input:  preprocessed/     msrmnts      [infield_steps]
    #                           basics       [ptrm_gained_vec]
    # output: preprocessed/     field_basics [Blab, field_dir_vec_initial, field_dir_vec]

    infield_steps = sc["preprocessed"]["msrmnts"]["infield_steps"] # works for Thermal-, Microwave- and Thellier- Thellier
    ptrm_gained_vec = sc["preprocessed"]["basics"]["ptrm_gained_vec"]

    # get lab field strength and direction of applied labfield
    Blab = infield_steps[0]["lab_field"]
    field_dir_vec_initial = helpers.dir2cart(infield_steps[0]["lab_field_dec"], infield_steps[0]["lab_field_inc"], 1)

    # check the field_dir_vec_initial, in an experiment the specimen could have been rotated so that field dir vec should be negative
    # check this by comparing the direction of the ptrm with field_dir_vec_initial
    vecl = ptrm_gained_vec[len(ptrm_gained_vec) - 1]    # last ptrm gained step

    vecl_abs = []                                       # make this vector absolute
    for i in range(len(vecl)):
        vecl_abs.append(abs(vecl[i]))

    # check field direction, and make it positive or negative is necessary, only works for fields applied along principal axis

    if field_dir_vec_initial[2] == 1 or field_dir_vec_initial[2] == -1:
        if field_dir_vec_initial[0] == 0 and field_dir_vec_initial[1] == 0:
            # then you have [0,0,1] or [0,0,-1]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [0, 0, -1]
            else:
                field_dir_vec = [0, 0, 1]
        else:       # [x,x,1] en x != 0
            field_dir_vec = field_dir_vec_initial
    elif field_dir_vec_initial[1] == 1 or field_dir_vec_initial[1] == -1:
        if field_dir_vec_initial[0] == 0 and field_dir_vec_initial[2] == 0:
            # then you have [0,1,0] or [0,-1,0]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [0, -1, 0]
            else:
                field_dir_vec = [0, 1, 0]
        else:       # [x,1,x] en x != 0
            field_dir_vec = field_dir_vec_initial
    elif field_dir_vec_initial[0] == 1 or field_dir_vec_initial[0] == -1:
        if field_dir_vec_initial[1] == 0 and field_dir_vec_initial[2] == 0:
            # then you have [1,0,0] or [-1,0,0]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [-1, 0, 0]
            else:
                field_dir_vec = [1, 0, 0]
        else:       # [1,x,x] with x != 0
            field_dir_vec = field_dir_vec_initial
    else:                                                                               # All other field directions that are not principle axis
        field_dir_vec = field_dir_vec_initial

    sc["preprocessed"]["field_basics"]["Blab"] = Blab
    sc["preprocessed"]["field_basics"]["field_dir_vec_initial"] = field_dir_vec_initial
    sc["preprocessed"]["field_basics"]["field_dir_vec"] = field_dir_vec
    return sc


def field_basics_pTh(sc):
    """
    Function to write the field basics to the suitcase for pseudo-Thellier. Only difference here w.r.t. the field_basics function above is that the "infield_steps" are now the ARM_aqc_steps for the NAA format, for pTh-GF (generic format) this is still the infieldstep. The lab field strength (Blab), and the direction of the labfield as provided by the user in the input (field_dir_vec_initial). The labfield direction is chcked with the direction of the last ptrm_step to check if the field_dir_vec_initial is correct or if it should be anti-parallel to that.

    input: the ptrm_gained_vec from the basisc, the infield_steps and the ARM_acq_steps from the measurements
    output: the field basics: Blab, field_dir_vec_initial, field_dir_vec
    """


    # input:  preprocessed/     msrmnts      [infield_steps, ARM_acq_steps]
    #                           basics       [ptrm_gained_vec]
    # output: preprocessed/     field_basics [Blab, field_dir_vec_initial, field_dir_vec]

    infield_steps = sc["preprocessed"]["msrmnts"]["infield_steps"]
    ARM_acq_steps = sc["preprocessed"]["msrmnts"]["ARM_acq_steps"]
    ptrm_gained_vec = sc["preprocessed"]["basics"]["ptrm_gained_vec"]

    # make extra if statement
    if len(infield_steps) == 0:
        infield_steps = ARM_acq_steps

    # get lab field strength and direction of applied labfield
    Blab = infield_steps[0]["lab_field"]
    field_dir_vec_initial = helpers.dir2cart(infield_steps[0]["lab_field_dec"], infield_steps[0]["lab_field_inc"], 1)

    # check the field_dir_vec_initial, in an experiment the specimen could have been rotated so that field dir vec should be negative
    # check this by comparing the direction of the ptrm with field_dir_vec_initial
    vecl = ptrm_gained_vec[len(ptrm_gained_vec) - 1]    # last ptrm gained step
    vecl_abs = []                                       # make this vector absolute
    for i in range(len(vecl)):
        vecl_abs.append(abs(vecl[i]))

    # check field direction, and make it positive or negative is necessary, only works for fields applied along principal axis
    if field_dir_vec_initial[2] == 1 or field_dir_vec_initial[2] == -1:
        if field_dir_vec_initial[0] == 0 and field_dir_vec_initial[1] == 0:
            # then you have [0,0,1] or [0,0,-1]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [0, 0, -1]
            else:
                field_dir_vec = [0, 0, 1]
        else:       # [x,x,1] en x != 0
            field_dir_vec = field_dir_vec_initial
    elif field_dir_vec_initial[1] == 1 or field_dir_vec_initial[1] == -1:
        if field_dir_vec_initial[0] == 0 and field_dir_vec_initial[2] == 0:
            # then you have [0,1,0] or [0,-1,0]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [0, -1, 0]
            else:
                field_dir_vec = [0, 1, 0]
        else:       # [x,1,x] en x != 0
            field_dir_vec = field_dir_vec_initial
    elif field_dir_vec_initial[0] == 1 or field_dir_vec_initial[0] == -1:
        if field_dir_vec_initial[1] == 0 and field_dir_vec_initial[2] == 0:
            # then you have [1,0,0] or [-1,0,0]
            if (max(vecl) - max(vecl_abs)) < 0:
                field_dir_vec = [-1, 0, 0]
            else:
                field_dir_vec = [1, 0, 0]
        else:       # [1,x,x] with x != 0
            field_dir_vec = field_dir_vec_initial
    else:                                                    # All other field directions that are not principal axis
        field_dir_vec = field_dir_vec_initial


    sc["preprocessed"]["field_basics"]["Blab"] = Blab
    sc["preprocessed"]["field_basics"]["field_dir_vec_initial"] = field_dir_vec_initial
    sc["preprocessed"]["field_basics"]["field_dir_vec"] = field_dir_vec
    return sc


def checks(sc):
    """
    Function to calculate the ptrm, tail, additivity checks, the vectors: [step,x,y,z], but also list of the moments for x_ptrm_check. The function is slit in two parts for the calculation of the modified thellier (thermal and microwave), and the Thellier-Thellier.

    input: specimen and specimen selection info, to know which steps are measured prior the check step, the measurement check data, the ptrm and NRM_rem vectors to calculate the gain or loss of a check
    output: the vectors and list for the different checks
    """


    # input:    input[specimen_selection, specimen]
    #           preprocessed/     basics[ptrm, NRM_rem]
    #           preprocessed/     msrmnts[ptrm_check_steps, tail_check_steps, add_check_steps, infield_antiparallel_steps]
    # output:   preprocessed/     checks[y_tail_check, x_tail_check, y_temp_tail_check, x_temp_ptrm_check, x_ptrm_check, y_ptrm_check, ptrm_check, add_check_step, add_check_vec, tail_check_vec, SCAT_ptrm_check_step, SCAT_tail_check_step]


    sam_sel = sc["input"]["specimen_selection"] # need this for ptrm check calculations, the order of steps is important here
    specimen = sc["input"]["specimen"] # need this for ptrm check calculations, the order of steps is important here


    ptrm_check_steps = sc["preprocessed"]["msrmnts"]["ptrm_check_steps"]
    tail_check_steps = sc["preprocessed"]["msrmnts"]["tail_check_steps"]
    add_check_steps = sc["preprocessed"]["msrmnts"]["add_check_steps"]
    infield_antiparallel_steps = sc["preprocessed"]["msrmnts"]["infield_antiparallel_steps"]

    ptrm = sc["preprocessed"]["basics"]["ptrm"]
    NRM_rem = sc["preprocessed"]["basics"]["NRM_rem"]


    y_tail_check = []
    y_temp_tail_check = []
    x_tail_check = []
    x_temp_ptrm_check = []
    x_ptrm_check = []
    y_ptrm_check = []
    SCAT_ptrm_check_step = []
    SCAT_tail_check_step = []

    tail_check_vec = []

    add_check_vec = []
    add_check_step = []

    ptrm_check_list= []

    def calc_mom(step, prevstep):
            moment = math.sqrt((step["x"] - prevstep["x"])**2 + (step["y"] - prevstep["y"])**2 + (step["z"] - prevstep["z"])**2)
            return moment

    def calc_mom_TT(step, prevstep):
            moment = helpers.norm(helpers.list_mult_num([ step["x"]-prevstep["x"], step["y"]- prevstep["y"], step["z"]-prevstep["z"]],0.5))
            return moment


    # "modified thellier" so Thermal Thellier of Microwave Thellier, is not calculated for pseudo Thellier
    if len(infield_antiparallel_steps) == 0:

        # Tail-check calculations
        for Tstep in tail_check_steps:
            y_tail_check.append(math.sqrt(Tstep["x"]**2 + Tstep["y"]**2 + Tstep["z"]**2))
            y_temp_tail_check.append(Tstep["step"])
            tail_check_vec.append([Tstep["x"],Tstep["y"],Tstep["z"]])

        # determine the x-coordinate for the Tail-check (need this for SCAT param)
        for check_step in y_temp_tail_check:
            for line in ptrm:                       # look through the calcualated ptrm gained specimen data
                if (line["step"] == check_step):    # find the corresponding ptrm-gained step with same "step" as the the check_step (this is temperature for Thermall Thellier)
                    x_tail_check.append(math.sqrt(line["x"]**2 + line["y"]**2 + line["z"]**2))  # calculate total magnetization and append
                    SCAT_tail_check_step.append(line["step"])

        # pTRM-check calculations
        for j in range(len(sam_sel)):
            if sam_sel[j]["type"] == 2:
                step = sam_sel[j]["step"]
                prevstep = sam_sel[j]["previous_step"]

                # Now find this step in the matrix with all the data to calculate x_ptrm_check & ptrm_check_list
                for i in range(len(specimen)):
                    if (specimen[i]["type"] == 2) and (specimen[i]["step"] == step) and (specimen[i]["previous_step"]==prevstep):
                        x_temp_ptrm_check.append(step)                            # copy treatment temp
                        x_ptrm_check.append(calc_mom(specimen[i], specimen[i - 1]))   # calculate the absolute gain in all three axis due to the pTRM-check step
                        ptrm_check_list.append([specimen[i]["step"], specimen[i]["x"] - specimen[i-1]["x"], specimen[i]["y"] - specimen[i-1]["y"], specimen[i]["z"] - specimen[i-1]["z"]])

        # determine the y-coordinate for the pTRM-check (need this for SCAT param)
        for check_step in x_temp_ptrm_check:
            for line in sam_sel:                    # look through all selected specimen data
                if (line["step"] == check_step) and (line["type"] == 0):   # find total_m for the corresponding zerofield step with same "step" as the the check_step (this is temperature for Thermall Thellier)
                    y_ptrm_check.append(line["total_m"])      # append the magnetization
                    SCAT_ptrm_check_step.append(line["step"])

        # additivity-check calculations
        for j in range(len(sam_sel)):
            if sam_sel[j]["type"] == 4:
                step = sam_sel[j]["step"]
                prevstep = sam_sel[j]["previous_step"]

                for i in range(len(specimen)):
                    if (specimen[i]["type"] == 4) and (specimen[i]["step"] == step) and (specimen[i]["previous_step"]==prevstep):
                        if specimen[i-1]["type"] == 1:       # additivity check after an infield/trm step
                            add_check_step.append([specimen[i]["step"], specimen[i-1]["step"]])  # copy treatment step and previous treatment step [i , j]
                            add_check_vec.append([specimen[i-1]["x"] - specimen[i]["x"], specimen[i-1]["y"] - specimen[i]["y"], specimen[i-1]["z"] - specimen[i]["z"]])
                        elif specimen[i-1]["type"] == 2:       # additivity check after an ptrm-check step
                            # adapted the formulas from SPD matlab code
                            add_check_step.append([specimen[i]["step"], specimen[i-2]["step"]])  # copy treatment step and previous treatment step [i , j]
                            add_check_vec.append([specimen[i-2]["x"] - specimen[i]["x"], specimen[i-2]["y"] - specimen[i]["y"], specimen[i-2]["z"] - specimen[i]["z"]])
                        elif (specimen[i-1]["type"] == 4) and (specimen[i-2]["type"] == 1) :       # additivity check after an aditivity check after an infield step
                            add_check_step.append([specimen[i]["step"], specimen[i-2]["step"]])  # copy treatment step and previous treatment step [i , j]
                            add_check_vec.append([specimen[i-2]["x"] - specimen[i]["x"], specimen[i-2]["y"] - specimen[i]["y"], specimen[i-2]["z"] - specimen[i]["z"]])


    elif len(infield_antiparallel_steps) != 0: # Thellier - Thellier
        # pTRM-check calculations
        for j in range(len(sam_sel)):
            if sam_sel[j]["type"] == 2:
                step = sam_sel[j]["step"]
                prevstep = sam_sel[j]["previous_step"]

                # Now find this step in the matrix with all the data to calculate x_ptrm_check & ptrm_check_list
                for i in range(len(specimen)):
                    if (specimen[i]["type"] == 2) and (specimen[i]["step"] == step) and (specimen[i]["previous_step"]==prevstep):
                        x_temp_ptrm_check.append(step)                            # copy treatment temp
                        x_ptrm_check.append(calc_mom_TT(specimen[i], specimen[i - 1]))  # calculate the absolute gain in all three axis due to the pTRM-check step
                        ptrm_check_list.append([specimen[i]["step"], (specimen[i]["x"] - specimen[i-1]["x"])*0.5, (specimen[i]["y"] - specimen[i-1]["y"])*0.5, (specimen[i]["z"] - specimen[i-1]["z"])*0.5])

        # determine the y-coordinate for the pTRM-check (need this for SCAT param)
        for check_step in x_temp_ptrm_check:
            for line in NRM_rem:       # look through NRM remaining data and find corresponding step
                if (line["step"] == check_step):
                    y_ptrm_check.append(math.sqrt(line["x"]**2 + line["y"]**2 + line["z"]**2))      # append the magnetization
                    SCAT_ptrm_check_step.append(line["step"])


    ptrm_check = []

    # rewrite format for ptrm_check
    if len(ptrm_check_list) != 0:
        for msrmt in ptrm_check_list:
            ptrm_check.append({
                "step": msrmt[0],
                "x": msrmt[1],
                "y": msrmt[2],
                "z": msrmt[3]
            })


    sc["preprocessed"]["checks"]["y_tail_check"] = y_tail_check
    sc["preprocessed"]["checks"]["x_tail_check"] = x_tail_check
    sc["preprocessed"]["checks"]["y_temp_tail_check"] = y_temp_tail_check
    sc["preprocessed"]["checks"]["x_temp_ptrm_check"] = x_temp_ptrm_check
    sc["preprocessed"]["checks"]["x_ptrm_check"] = x_ptrm_check
    sc["preprocessed"]["checks"]["y_ptrm_check"] = y_ptrm_check
    sc["preprocessed"]["checks"]["ptrm_check"] = ptrm_check
    sc["preprocessed"]["checks"]["add_check_step"] = add_check_step
    sc["preprocessed"]["checks"]["add_check_vec"] = add_check_vec
    sc["preprocessed"]["checks"]["tail_check_vec"]= tail_check_vec
    sc["preprocessed"]["checks"]["SCAT_ptrm_check_step"]= SCAT_ptrm_check_step
    sc["preprocessed"]["checks"]["SCAT_tail_check_step"]= SCAT_tail_check_step
    return sc
