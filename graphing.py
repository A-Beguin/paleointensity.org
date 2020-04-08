import math
import helpers

# functions to make data for the four plots:
# plot_arai, plot_zijderveld, plot_equalarea, plot_magnitisation
# all functions are suitable for the classical Thellier-Thellier and for modified Thellier, both for Thermal and Microwave

# returns [ptrm_gained, nrm_rem, ptrmCheck, tailCheck, order_steps]
# ptrm      = [tmp, x], x value per tmp step
# nrm       = [tmp, y], y value per tmp step
# ptrmCheck = [tmp, x, y, xprevious, yprevious]
# tailCheck = [tmp, x, y]


def plot_arai(data):
    ptrm_gained = []
    nrm_rem = []
    ptrmCheck = []
    tailCheck = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))
    infield_antiparallel_steps = list(filter(lambda msrmt: msrmt['type'] == 5, data))


    NRMrem_list = []
    ptrm_list = []


    # check which format you have, "classical" or "modified" Thellier
    if len(infield_antiparallel_steps) != 0:    # "classical thellier Thelier-Thellier (TT)"
        # first the NRM, which is always a zerofield step, also for the Thellier - Thellier protocol, might be that this is not given so check if present
        if len(zerofield_steps) != 0:
            NRMrem_list.append([zerofield_steps[0]["step"], zerofield_steps[0]["x"], zerofield_steps[0]["y"], zerofield_steps[0]["z"]])
            ptrm_list.append([zerofield_steps[0]["step"], 0, 0, 0])  # dit zorgt voor je punt op de y as.

        # determine the NRM remaining & ptrm gained
        for i_step in infield_steps:
            for ai_step in infield_antiparallel_steps:
                if i_step["step"] == ai_step["step"]:
                    NRMrem_list.append([i_step["step"], (i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.])
                    ptrm_list.append([i_step["step"], (i_step["x"] - ai_step["x"]) / 2., (i_step["y"] - ai_step["y"]) / 2., (i_step["z"] - ai_step["z"]) / 2.])

    if len(infield_antiparallel_steps) == 0: # "modified thellier"

        for i in range(1):
            ptrm_list.append([zerofield_steps[i]["step"], 0, 0, 0])    # y-axis first point.
            NRMrem_list.append([zerofield_steps[i]["step"], zerofield_steps[i]["x"], zerofield_steps[i]["y"], zerofield_steps[i]["z"]]) # added this line 05-02

        # make ptrm_list & NRM_list when the steps are the same
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])
                    ## make NRMrem_list
                    #for z_step in zerofield_steps: (05-02 changed!!)
                    NRMrem_list.append([z_step["step"], z_step["x"], z_step["y"], z_step["z"]])

    # get NRM0
    nrm0 = math.sqrt(NRMrem_list[0][1]**2 + NRMrem_list[0][2]**2 + NRMrem_list[0][3]**2)

    # rewrite format ptrm and NRM
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

    nrm_rem = []
    for z_step in NRM_rem:
        nrm_rem.append([z_step["step"], (math.sqrt(z_step["x"]**2 + z_step["y"]**2 + z_step["z"]**2)) / nrm0])

    ptrm_gained = []
    for step in ptrm:
        ptrm_gained.append([step["step"], (math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2)) / nrm0])

    # helper function
    def calc_mom(step, laststep):
        moment = math.sqrt((step["x"] - laststep["x"])**2 + (step["y"] - laststep["y"])**2 + (step["z"] - laststep["z"])**2)
        return moment

    def list_mult_num(l, num):
        result = []
        for i in range(len(l)):
            result.append(l[i] * num)
        return result

    def calc_mom_TT(step, prevstep):
        moment = helpers.norm(list_mult_num([ step["x"]-prevstep["x"], step["y"]- prevstep["y"], step["z"]-prevstep["z"]],0.5))
        return moment


    # pTRM-check calculations
    temppTRM = []
    x = []
    y = []
    xprev = []
    yprev = []
    temp_prev = []

    if len(infield_antiparallel_steps) == 0: # Modified Thellier protocol use calc_mom
        for j in range(len(data)):
            if data[j]["type"] == 2:
                temppTRM.append(data[j]["step"])  # copy treatment temp
                x.append(calc_mom(data[j], data[j - 1])/ nrm0)  # calculate the absolute gain in all three axis due to the pTRM-check step
                temp_prev.append(data[j]["previous_step"])

    if len(infield_antiparallel_steps) != 0: # Thellier - Thellier use calc_mom_TT
        for j in range(len(data)):
            if data[j]["type"] == 2:
                temppTRM.append(data[j]["step"])                   # copy treatment temp
                x.append(calc_mom_TT(data[j], data[j - 1])/ nrm0)  # calculate the absolute gain in all three axis due to the pTRM-check step
                temp_prev.append(data[j]["previous_step"])

    for Tstep in temppTRM:
        for z_step in NRM_rem:
            if z_step["step"] == Tstep:
                y.append((math.sqrt(z_step["x"]**2 + z_step["y"]**2 + z_step["z"]**2)) / nrm0)      # find the NRM rem of the corresponding step

    for preTstep in temp_prev:
        for z_step in NRM_rem:
            if z_step["step"] == preTstep:
                yprev.append((math.sqrt(z_step["x"]**2 + z_step["y"]**2 + z_step["z"]**2)) / nrm0)  # find the NRM rem of the corresponding previous step

    for preTstep in temp_prev:
        for step in ptrm:
            if step["step"] == preTstep:
                xprev.append(math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2) / nrm0)  # find the ptrm gained of the corresponding previous step


    # print("x",x)
    # print("y",y)
    # print("xprev",xprev)
    # print("yprev",yprev)
    for le in range(len(temppTRM)):
        ptrmCheck.append([temppTRM[le], x[le], y[le], xprev[le], yprev[le]])  # the above calculated values to the ptrmCheck array


    # Tail-check calculations : tailCheck = [tmp, x, y]
    y_tail_check = []
    step_tailCheck = []

    for j in range(len(data)):
        if data[j]["type"] == 3:                       # Then it is a tail check!
            step_tailCheck.append(data[j]["step"])          # copy treatment step
            y_tail_check.append(data[j]["total_m"] / nrm0)  # Append the y coord of the tail check, this is magnetization scaled by NRM0

    # determine the x-coordinate for the Tail-check
    x_tail_check = []
    for check_step in step_tailCheck:       # go though steps ("temperature", or "AGV integral")
        for line in ptrm:                   # look trough the calculated ptrm gained specimen data
            if (line["step"] == check_step):  # find the corresponding ptrm-gained step with same "step" as the the check_step (this is temperature for Thermall Thellier)
                x_tail_check.append(math.sqrt(line["x"]**2 + line["y"]**2 + line["z"]**2) / nrm0)  # calculate total magnetization, scale to NRM0, and append

    for le in range(len(step_tailCheck)):
        tailCheck.append([step_tailCheck[le], x_tail_check[le], y_tail_check[le]])  # the above calutated vaules to the Tail-Check array


    # determine IZ or ZI order
    # for TT format IZ means first parallel infield, then antiparrallel, ZI means first aniparallel then paralell.
    order_steps = []
    ind_NRM  = []
    for nstep in NRM_rem:
        for i in range(len(data)):
            if (nstep["step"] == data[i]["step"]) and (data[i]["type"] == 0 or data[i]["type"] == 5):
                ind_NRM.append([nstep["step"], i])

    ind_pTRM  = []
    for pstep in ptrm:
        for i in range(len(data)):
            if (pstep["step"] == data[i]["step"]) and (data[i]["type"] == 1):
                ind_pTRM.append([pstep["step"], i])

    for nstep in ind_NRM:
        for pstep in ind_pTRM:
            if (nstep[0] == pstep[0]):      # same temp check index
                if (nstep[1] < pstep[1]):
                    order_steps.append([nstep[0], "ZI"]) # index of the NRM smaler than the ptrm then first zero then infield
                elif (nstep[1] > pstep[1]):
                    order_steps.append([nstep[0], "IZ"])

    return[ptrm_gained, nrm_rem, ptrmCheck, tailCheck, order_steps]


#returns [horizontal_NorthUP, vertical_NorthUP, horizontal_WestUP, vertical_WestUP]
# horizontal = [temp, x, y], the horizontal projection
# vertical   = [temp, x, y], the vertical projection
def plot_zijderveld(data):
    horizontal_NorthUP = []
    vertical_NorthUP = []
    horizontal_WestUP = []
    vertical_WestUP = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))      # only one needed for "modified Thellier"

    # for Classical Thellier need both infield and antiparrallel infield steps
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))
    infield_antiparallel_steps = list(filter(lambda msrmt: msrmt['type'] == 5, data))


    NRMrem_list = []
    # check which format you have, "classical" or "modified" Thellier
    if len(infield_antiparallel_steps) != 0:    # "classical thellier Thelier-Thellier (TT)"
        # first the NRM, which is always a zerofield step, also for the Thellier - Thellier protocol, migth be that this is not given so check if present
        if len(zerofield_steps) != 0:
            NRMrem_list.append([zerofield_steps[0]["step"], zerofield_steps[0]["x"], zerofield_steps[0]["y"], zerofield_steps[0]["z"]])

        # determine the NRM remaining & ptrm gained for the selection
        for i_step in infield_steps:
            for ai_step in infield_antiparallel_steps:
                if i_step["step"] == ai_step["step"]:
                    NRMrem_list.append([i_step["step"], (i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.])

    if len(infield_antiparallel_steps) == 0: # "modified thellier"
        # make NRMrem_list
        for z_step in zerofield_steps:
            NRMrem_list.append([z_step["step"], z_step["x"], z_step["y"], z_step["z"]])

    NRM_rem = []
    for msrmt in NRMrem_list:
        NRM_rem.append({
            "step": msrmt[0],
            "x": msrmt[1],
            "y": msrmt[2],
            "z": msrmt[3]
        })

    for z_step in NRM_rem:
        horizontal_NorthUP.append([z_step["step"], z_step["y"], z_step["x"]])
        vertical_NorthUP.append([z_step["step"], z_step["y"], -1 * z_step["z"]])
        horizontal_WestUP.append([z_step["step"], z_step["x"], -1 * z_step["y"]])
        vertical_WestUP.append([z_step["step"], z_step["x"], -1 * z_step["z"]])

    return [horizontal_NorthUP, vertical_NorthUP, horizontal_WestUP, vertical_WestUP]


#returns [deMag, gainedMag]
# deMag     = [temp, dec, inc], direction per temp zero-step
# gainedMag = [temp, dec, inc], direction per temp partial TRM gained step
def plot_equalarea(data):
    deMag = []
    gainedMag = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))

    infield_antiparallel_steps = list(filter(lambda msrmt: msrmt['type'] == 5, data)) # needed for TT

    # helper functions, calulate Dec & Inc
    def dec(step):
        declination = math.atan2(step["y"], step["x"]) * 180 / math.pi
        return declination

    def inc(step):
        inclination = math.asin(step["z"] / (math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2))) * 180 / math.pi
        return inclination

    # make Demay array and pTRM-gained array

    # first make NRM_rem & ptrm

    NRMrem_list = []
    ptrm_list = []
    # check which format you have, "classical" or "modified" Thellier
    if len(infield_antiparallel_steps) != 0:    # "classical thellier Thelier-Thellier (TT)"
        # first the NRM, which is always a zerofield step, also for the Thellier - Thellier protocol, might be that this is not given so check if present
        if len(zerofield_steps) != 0:
            NRMrem_list.append([zerofield_steps[0]["step"], zerofield_steps[0]["x"], zerofield_steps[0]["y"], zerofield_steps[0]["z"]])

        # determine the NRM remaining & ptrm gained
        for i_step in infield_steps:
            for ai_step in infield_antiparallel_steps:
                if i_step["step"] == ai_step["step"]:
                    NRMrem_list.append([i_step["step"], (i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.])
                    ptrm_list.append([i_step["step"], (i_step["x"] - ai_step["x"]) / 2., (i_step["y"] - ai_step["y"]) / 2., (i_step["z"] - ai_step["z"]) / 2.])

    if len(infield_antiparallel_steps) == 0: # "modified thellier"
        # make NRMrem_list
        for z_step in zerofield_steps:
            NRMrem_list.append([z_step["step"], z_step["x"], z_step["y"], z_step["z"]])

        # for the gainedMag first calcuate the pTRM-gained
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

    # rewrite format ptrm and NRM_rem
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


    # fill the deMag array
    for z_step in NRM_rem:
        deMag.append([z_step["step"], dec(z_step), inc(z_step)])

    # fill the gainedMag array with same functions for dec & inc
    for ptrm_step in ptrm:
        gainedMag.append([ptrm_step["step"], dec(ptrm_step), inc(ptrm_step)])

    return [deMag, gainedMag]


#returns [m_deMag, m_gainedMag]
# m_deMag      = [temp, m], temp=x, m=y, moment for zero-field steps
# m_gainedMag   = [temp, m], temp=x, m=y, moment for partial TRM gained
def plot_magnitisation(data):
    m_deMag = []
    m_gainedMag = []

    zerofield_steps = list(filter(lambda msrmt: msrmt['type'] == 0, data))
    infield_steps = list(filter(lambda msrmt: msrmt['type'] == 1, data))

    infield_antiparallel_steps = list(filter(lambda msrmt: msrmt['type'] == 5, data))


    def calc_m(step):
        moment = math.sqrt(step["x"]**2 + step["y"]**2 + step["z"]**2)
        return moment


    # start COPY from ARAI PLOT function
    NRMrem_list = []
    ptrm_list = []

    # check which format you have, "classical" or "modified" Thellier
    if len(infield_antiparallel_steps) != 0:    # "classical thellier Thelier-Thellier (TT)"
        # first the NRM, which is always a zerofield step, also for the Thellier - Thellier protocol, might be that this is not given so check if present
        if len(zerofield_steps) != 0:
            NRMrem_list.append([zerofield_steps[0]["step"], zerofield_steps[0]["x"], zerofield_steps[0]["y"], zerofield_steps[0]["z"]])

        # determine the NRM remaining & ptrm gained for the selection
        for i_step in infield_steps:
            for ai_step in infield_antiparallel_steps:
                if i_step["step"] == ai_step["step"]:
                    NRMrem_list.append([i_step["step"], (i_step["x"] + ai_step["x"]) / 2., (i_step["y"] + ai_step["y"]) / 2., (i_step["z"] + ai_step["z"]) / 2.])
                    ptrm_list.append([i_step["step"], (i_step["x"] - ai_step["x"]) / 2., (i_step["y"] - ai_step["y"]) / 2., (i_step["z"] - ai_step["z"]) / 2.])

    if len(infield_antiparallel_steps) == 0: # "modified thellier"
        # make ptrm_list
        for i_step in infield_steps:
            for z_step in zerofield_steps:
                if i_step["step"] == z_step["step"]:
                    ptrm_list.append([i_step["step"], (i_step["x"] - z_step["x"]), (i_step["y"] - z_step["y"]), (i_step["z"] - z_step["z"])])

        # make NRMrem_list
        for z_step in zerofield_steps:
            NRMrem_list.append([z_step["step"], z_step["x"], z_step["y"], z_step["z"]])

    # get NRM0
    nrm0 = math.sqrt(NRMrem_list[0][1]**2 + NRMrem_list[0][2]**2 + NRMrem_list[0][3]**2)

    # rewrite format ptrm and NRM
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
    # stop COPY from ARAI PLOT function


    m_deMag = []
    for z_step in NRM_rem:
        m_deMag.append([z_step["step"], (math.sqrt(z_step["x"]**2 + z_step["y"]**2 + z_step["z"]**2)) / nrm0])


    m_gainedMag = [[NRM_rem[0]["step"], 0]]
    for p_step in ptrm:
        m_gainedMag.append([p_step["step"], (math.sqrt(p_step["x"]**2 + p_step["y"]**2 + p_step["z"]**2)) / nrm0])

    return [m_deMag, m_gainedMag]
