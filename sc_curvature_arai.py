import numpy
import math
import helpers_curv


def AraiCurvature(sc):
    """
    Function for calculating the radius of the best fit circle to a set of
    x-y coordinates when there are at least three points available. This uses the curvature helpers function (helpers_curv) in which all the calculations are done.

    input from suitcase: list of x points, list of y points
    output to suitcase: k, SSE.  curvature and SSE

    """
    # input:    preprocessed/           basics[x_ptrm, y_nrm]
    # output:   curv_arai_statistics/   AraiCurvature[k, SSE]

    x = sc["preprocessed"]["basics"]["x_ptrm"]
    y = sc["preprocessed"]["basics"]["y_nrm"]

    if len(x)>3:
        (k, SSE) = helpers_curv.kprime_calc(x,y)

        sc["curv_arai_statistics"]["AraiCurvature"]["k_prime"] = k
        sc["curv_arai_statistics"]["AraiCurvature"]["SSE_k_prime"] = SSE
    return sc
