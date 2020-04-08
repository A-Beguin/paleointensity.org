import math

# this is partly a copy of the graphing.py document. Change the work flow that it is suitable for both pTh-GF & pTh-NAA

# returns [ptrm_gained, nrm_rem, order_steps]
# ptrm      = [tmp, x], x value per tmp step
# nrm       = [tmp, y], y value per tmp step

# optimized for the use of both GF and NAA for pTh


def plot_arai_pTh(data):
    ptrm_gained = []
    nrm_rem = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))    # for GF

    ARM_acq_steps_all = list(filter(lambda m: m['type'] == 6, data))        # for NAA

    # first calculate the pTRM-gained for x, y and z individually
    # assumption here that the extra zero field steps are always before the start of the infield steps!
    ptrm_list = []

    if len(infield_steps) != 0:  #pseudoThellier - GF
        for i in range(len(zerofield_steps) - len(infield_steps)):
            ptrm_list.append([zerofield_steps[i]["step"], 0, 0, 0])
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

    if len(ARM_acq_steps_all) != 0:  #pseudoThellier - NAA , Correct for the NRM remaining which is the first ARM aqc step
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
        for a_step in ARM_acq_steps_all:
            ptrm_list.append([a_step["step"], a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0])

    # same format
    ptrm = []
    for msrmt in ptrm_list:
        ptrm.append({
            "step": msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })

    nrm0 = zerofield_steps[0]["total_m"]

    nrm_rem = []
    for z_step in zerofield_steps:
        nrm_rem.append([z_step["step"], z_step["total_m"] / nrm0])

    ptrm_gained = []
    for step in ptrm:
        ptrm_gained.append([step["step"], (math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2)) / nrm0])

    # determine IZ or ZI order
    order_steps = []
    for step in zerofield_steps:
        if float(step["previous_step"]) == float(step["step"]):
            order_steps.append([step["step"], "IZ"])
        elif float(step["previous_step"]) < float(step["step"]):
            order_steps.append([step["step"], "ZI"])

    # empty list are for website application (checks lists )
    return[ptrm_gained, nrm_rem, [], [], order_steps]


#returns [horizontal_NorthUP, vertical_NorthUP, horizontal_WestUP, vertical_WestUP]
# horizontal = [temp, x, y], horizontal projection
# vertical   = [temp, x, y], vertical projection

# nothing changed for the pTh and other (Thermal)Thellier graphing functions
def plot_zijderveld_pTh(data):
    horizontal_NorthUP = []
    vertical_NorthUP = []
    horizontal_WestUP = []
    vertical_WestUP = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))

    for z_step in zerofield_steps:
        horizontal_NorthUP.append([z_step["step"], z_step["y"], z_step["x"]])
        vertical_NorthUP.append([z_step["step"], z_step["y"], -1 * z_step["z"]])
        horizontal_WestUP.append([z_step["step"], z_step["x"], -1 * z_step["y"]])
        vertical_WestUP.append([z_step["step"], z_step["x"], -1 * z_step["z"]])

    return [horizontal_NorthUP, vertical_NorthUP, horizontal_WestUP, vertical_WestUP]


#returns [deMag, gainedMag]
# deMag     = [temp, dec, inc], direction per temp zero-step
# gainedMag = [temp, dec, inc], direction per temp partial TRM gained step

# optimized for the use of both GF and NAA for pTh
def plot_equalarea_pTh(data):
    deMag = []
    gainedMag = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))    # for pTh

    ARM_acq_steps_all = list(filter(lambda m: m['type'] == 6, data))        # for NAA

    # helper functions, calculate Dec & Inc
    def dec(step):
        declination = math.atan2(step["y"], step["x"]) * 180 / math.pi
        return declination

    def inc(step):
        inclination = []
        R = math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2)
        if R > 0:
            inclination = math.asin(step["z"] / (R)) * 180 / math.pi
        return inclination

    # fill the deMag array, this one is easy, just the zerofieldsteps
    for z_step in zerofield_steps:
        deMag.append([z_step["step"], dec(z_step), inc(z_step)])

    # for the gainedMag first calculate the pTRM-gained
    # for pTh_NAA the ptrm_gained (Arm acquisition) is calculated differently, check if you have GF or NAA
    ptrm_list = []
    if len(infield_steps) == 0:             # then NAA
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
        for a_step in ARM_acq_steps_all:
            if int(a_step["step"]) > 0:          # first step has no field
                ptrm_list.append([a_step["step"], a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0])

    else:                                   # pTh-GF
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

    # same format as data
    ptrm = []
    for msrmt in ptrm_list:
        ptrm.append({
            "step": msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })

    # fill the gainedMag array with functions for dec & inc
    for ptrm_step in ptrm:
        gainedMag.append([ptrm_step["step"], dec(ptrm_step), inc(ptrm_step)])

    return [deMag, gainedMag]


#returns [m_deMag, m_gainedMag]
# m_deMag      = [temp, m], temp=x, m=y, moment for zero-field steps
# m_gainedMag   = [temp, m], temp=x, m=y, moment for partial TRM gained

# optimized for the use of both GF and NAA for pTh
def plot_magnitisation_pTh(data):
    m_deMag = []
    m_gainedMag = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))    # for pTh

    ARM_acq_steps_all = list(filter(lambda m: m['type'] == 6, data))        # for NAA

    def calc_m(step):
        moment = math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2)
        return moment

    NRM0 = math.sqrt(zerofield_steps[0]["x"]**2 + zerofield_steps[0]["y"]**2 + zerofield_steps[0]["z"]**2)

    # fill the deMag array, this one is easy, just the zerofieldsteps
    for z_step in zerofield_steps:
        m_deMag.append([z_step["step"], calc_m(z_step) / NRM0])

    # first vaule for the gained array should be at the same temp as NRM and moment = 0
    m_gainedMag = [[zerofield_steps[0]["step"], 0]]

    # to fill rest of the gained array first check if you have GF or NAA
    ptrm_list = []
    if len(infield_steps) == 0:  # then NAA
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
         # fill the m_gainedMag array
        for a_step in ARM_acq_steps_all:
            ptrm_list.append([a_step["step"], a_step["x"] - x0, a_step["y"] -y0 , a_step["z"] - z0])
    else:  # then GF
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

    # same format as data
    ptrm = []
    for msrmt in ptrm_list:
        ptrm.append({
            "step": msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })
    # stop COPY from equal area function

    # fill the m_gainedMag array with same functions for calculation of moment
    for ptrm_step in ptrm:
        m_gainedMag.append([ptrm_step["step"], calc_m(ptrm_step) / NRM0])

    return [m_deMag, m_gainedMag]


#returns [ARM_ARM, Dem_Dem]
# Dem_Dem       = [step, NRMdem, ARMdem]  #=(step, x,y)
# ARM_ARM       = [step, ARMacq, ARMdem]  #=(step, x,y)

def plot_pTh_graphs(data):
    ARM_ARM = []
    Dem_Dem = []

    NRM_dem_steps = list(filter(lambda m: m['type'] == 0, data))
    ARM_acq_steps_all = list(filter(lambda m: m['type'] == 6, data))
    ARM_dem_steps = list(filter(lambda m: m['type'] == 7, data))

    nrm0 = NRM_dem_steps[0]["total_m"]


    if len(ARM_acq_steps_all) != 0:
        arm = []
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
        for a_step in ARM_acq_steps_all:
            arm.append([a_step["step"], math.sqrt((a_step["x"] - x0)**2 + (a_step["y"] -y0)**2 + (a_step["z"] - z0)**2)])

    if len(ARM_dem_steps) != 0:
        arm_dem = []
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
        for ad_step in ARM_dem_steps:
            arm_dem.append([ad_step["step"], math.sqrt((ad_step["x"] - x0)**2 + (ad_step["y"] -y0)**2 + (ad_step["z"] - z0)**2)])

    # only if the ARM demag is given try to make data for the plots, otherwise return empty lists
    if len(ARM_dem_steps) != 0:
        for n_step in NRM_dem_steps:
            for i in range(len(arm_dem)):
                if n_step["step"] == arm[i][0]:
                    Dem_Dem.append([n_step["step"], n_step["total_m"] / nrm0, arm_dem[i][1] / nrm0])
        for i in range(len(arm)):
            for d_step in ARM_dem_steps:
                if arm[i][0] == d_step["step"]:
                    ARM_ARM.append([d_step["step"], arm[i][1] / nrm0, d_step["total_m"] / nrm0])

    return [ARM_ARM, Dem_Dem]


#returns [B12_ARM]
# B12_ARM      = [B12_ARM]

# function to calculate the B12ARM value using all ARM acquisition steps, this can be separate from the suitcase since this parameters only needs to be calculated once.
def B12ARM_param_pth(data):

    ARM_acq_steps_all = list(filter(lambda m: m['type'] == 6, data))

    if len(ARM_acq_steps_all) == 0:  # then you probably have pTh-GF # copy code from Arai-plot graphing
        zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
        infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))

        # first calculate the pTRM-gained for x, y and z individually
        # assumption here that the extra zero field steps are always before the start of the infield steps!
        ptrm_list = []
        for i in range(len(zerofield_steps) - len(infield_steps)):
            ptrm_list.append([zerofield_steps[i]["step"], 0, 0, 0])
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

        # same format
        ptrm = []
        for msrmt in ptrm_list:
            ptrm.append({
                "step": msrmt[0],
                "x": msrmt[1],
                "y": msrmt[2],
                "z": msrmt[3]
            })

        ptrm_gained = []
        for step in ptrm:
            ptrm_gained.append([step["step"], (math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2))])

        # end copy data

        # Make two list with the ARM, and the Field for all steps
        ARM_list = []
        Field_list = []
        for i in range(len(ptrm_gained)):
            ARM_list.append(ptrm_gained[i][1])
            Field_list.append(ptrm_gained[i][0])

    else:                           # so not zero, then you have the pTh-NAA format and calculations are easy
        # Make two list with the ARM, and the Field for all steps
        ARM_list = []
        Field_list = []
        arm = []
        x0 = ARM_acq_steps_all[0]["x"]
        y0 = ARM_acq_steps_all[0]["y"]
        z0 = ARM_acq_steps_all[0]["z"]
        for a_step in ARM_acq_steps_all:
            ARM_list.append(math.sqrt((a_step["x"] - x0)**2 + (a_step["y"] -y0)**2 + (a_step["z"] - z0)**2))
            Field_list.append(a_step["step"])

    # now we have the ARM and Field list, works for both GF and NAA format.

    # saturation ARM is last value in the ARM acq list
    SARM = ARM_list[len(ARM_list) - 1]

    # find halve of that value
    SARM_12 = 0.5 * SARM

    # determine the points just above and below the ARMsat
    # pay attention to the upper and lower bounds, if there is a case that [ind_min+1] does not exist, do not calculate B12ARM

    ind_min = 999
    for i in range(len(ARM_list) - 1):
        if (ARM_list[i] <= SARM_12) & (ARM_list[i + 1] > SARM_12):
            ind_min = i

    if ind_min == 999:
        B12_ARM = None
    else:
        Alow = ARM_list[ind_min]
        Aup = ARM_list[ind_min + 1]

        Flow = Field_list[ind_min]
        Fup = Field_list[ind_min + 1]

        # calculate the line formula between these two points
        a = (Aup - Alow) / (Fup - Flow)
        b = Alow - a * Flow

        # determine the field value (x) for the known halve Saturation ARM
        B12_ARM = (SARM_12 - b) / a

    return [B12_ARM]
