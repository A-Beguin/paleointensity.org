import math
import numpy
import helpers_curv
import helpers


"""
This script contains the helper functions for the calculation of the Arai-plot curvature.
"""

def kprime_calc(x, y):
    """
    This function is the main body of the k_prime calculation and uses the other helper functions in this script. This function calculates the radius of the best fit circle to a set of
    x-y coordinates.

    Paterson, G. A., (2011), A simple test for the presence of multidomain
    behaviour during paleointensity experiments, J. Geophys. Res., in press,
    doi: 10.1029/2011JB008369


    input: list of x points, list of y points
    output: k, SSE.  curvature and SSE
    """

    # Normalize the x and y points by the maximum
    n = len(x)
    X = []
    Y = []
    for i in range(n):
        X.append(x[i] / max(x))
        Y.append(y[i] / max(y))

    # Determine the initial estimate
    E1 = helpers_curv.TaubinSVD(X, Y)

    a = E1[0]
    b = E1[1]
    r = E1[2]

    if n > 3:

        E2 = helpers_curv.LMA(X, Y, E1)

        a = E2[0]
        b = E2[1]
        r = E2[2]

        # Determine the sense of curvature
        if a <= sum(X) / len(X) and b <= sum(Y) / len(Y):
            k = -1. / r
        else:
            k = 1. / r

        SSE = helpers_curv.get_SSE(a, b, r, X, Y)

    return(k, SSE)



def TaubinSVD(X, Y):
    """
    Algebraic circle fit by Taubin, helper function for calculating the curvature

    G. Taubin, "Estimation Of Planar Curves, Surfaces And Nonplanar
                  Space Curves Defined By Implicit Equations, With
                  Applications To Edge And Range Image Segmentation",
      IEEE Trans. PAMI, Vol. 13, pages 1115-1138, (1991)


    Input: lists of x values and y values
    output: a, b, r.  a and b are the center of the fitting circle, and r is the radius
    """

    # Get the number of points
    n = len(X)

    # Center the XY data
    centroid = [sum(X) / n, sum(Y) / n]

    Z = []
    tmp_X = []
    tmp_Y = []
    for i in range(n):
        tmp_X.append(X[i] - centroid[0])
        tmp_Y.append(Y[i] - centroid[1])
        Z.append(tmp_X[i] * tmp_X[i] + tmp_Y[i] * tmp_Y[i])

    Zmean = sum(Z) / n

    Z0 = []
    for i in range(n):
        Z0.append((Z[i] - Zmean) / (2. * math.sqrt(Zmean)))

    ZXY = [Z0, tmp_X, tmp_Y]

    # rewrtite with numpy

    ZXY = helpers.transpose_list(ZXY)

    # 1) numpy arrays
    ZXY = numpy.array(ZXY)

    # 2) do numpy, Do the SVD decomposition
    U, S, V = numpy.linalg.svd(ZXY, full_matrices=False)

    # 3) numpy back to python lists
    V = V.tolist()
    U = U.tolist()
    S = S.tolist()

    A = V[2]

    A[0] = A[0] / (2. * math.sqrt(Zmean))
    A.append(-1. * Zmean * A[0])

    # Get the circle center (a,b)
    a = (-1 * A[1]) / A[0] / 2 + centroid[0]
    b = (-1 * A[2]) / A[0] / 2 + centroid[1]

    # Circle radius
    # AB added abs () in sqrt and commented the ; at the end of the line
    r = math.sqrt(abs(A[1] * A[1] + A[2] * A[2] - 4 * A[0] * A[3])) / abs(A[0]) / 2  # ;

    return(a, b, r)


def VarCircle(X, Y, ParIn):
    """
    Function that computes the specimen variance of distances from data points (XY) to the circle Par = [a b R]
    helper function for calculating the curvature

    Input: list of of x and y values, and a tuple containing the parameters of the circle (a, b, r)

    """

    # Handle inputs
    n = len(X)

    Dx = []
    Dy = []
    D = []
    for i in range(n):
        Dx.append(X[i] - ParIn[0])
        Dy.append(Y[i] - ParIn[1])
        D.append(math.sqrt(Dx[i] * Dx[i] + Dy[i] * Dy[i]) - ParIn[2])

    result = helpers.dot_product(D, D) / (n - 3)

    return(result)


def LMA(X, Y, ParIni):
    """
    Geometric circle fit (minimizing orthogonal distances), helper function for calculating the curvature

    based on the Levenberg-Marquardt scheme in the "algebraic parameters"
    A,B,C,D  with constraint B*B+C*C-4*A*D=1

        N. Chernov and C. Lesort, "Least squares fitting of circles",
        J. Math. Imag. Vision, Vol. 23, 239-251 (2005)

    input: lists x values (X) and y values (Y), and a tuple containing an initial guess (a, b, r)
           which is acquired by using an algebraic circle fit (TaubinSVD)

    output: a, b, r.  a and b are the center of the fitting circle, and r is the radius
    """

    # Get the number of points
    n = len(X)

    # The algorithm parameters
    factorUp = 10
    factorDown = 0.04
    lambda0 = 0.01
    epsilon = 0.000001
    IterMAX = 50
    AdjustMax = 20
    Xshift = 0
    Yshift = 0
    dX = 1
    dY = 0

    anew = ParIni[0] + Xshift
    bnew = ParIni[1] + Yshift
    Anew = 1. / (2. * ParIni[2])
    aabb = anew * anew + bnew * bnew
    Fnew = (aabb - ParIni[2] * ParIni[2]) * Anew
    Tnew = math.acos(-anew / math.sqrt(aabb))

    if bnew > 0:
        Tnew = 2 * math.pi - Tnew
    VarNew = helpers_curv.VarCircle(X, Y, ParIni)

    VarLambda = lambda0
    finish = 0

    for it in range(0, IterMAX):

        Aold = Anew
        Fold = Fnew
        Told = Tnew
        VarOld = VarNew

        H = math.sqrt(1 + 4 * Aold * Fold)
        aold = -H * math.cos(Told) / (Aold + Aold) - Xshift
        bold = -H * math.sin(Told) / (Aold + Aold) - Yshift
        Rold = 1 / abs(Aold + Aold)

        DD = 1 + 4 * Aold * Fold
        D = math.sqrt(DD)
        CT = math.cos(Told)
        ST = math.sin(Told)

        # Initialize the sums
        H11 = 0
        H12 = 0
        H13 = 0
        H22 = 0
        H23 = 0
        H33 = 0
        F1 = 0
        F2 = 0
        F3 = 0

        for i in range(0, n):

            Xi = X[i] + Xshift
            Yi = Y[i] + Yshift
            Zi = Xi * Xi + Yi * Yi
            Ui = Xi * CT + Yi * ST
            Vi = -Xi * ST + Yi * CT

            ADF = Aold * Zi + D * Ui + Fold
            SQ = math.sqrt(4 * Aold * ADF + 1)
            DEN = SQ + 1
            Gi = 2 * ADF / DEN
            FACT = 2 / DEN * (1 - Aold * Gi / SQ)
            DGDAi = FACT * (Zi + 2 * Fold * Ui / D) - Gi * Gi / SQ
            DGDFi = FACT * (2 * Aold * Ui / D + 1)
            DGDTi = FACT * D * Vi

            H11 = H11 + DGDAi * DGDAi
            H12 = H12 + DGDAi * DGDFi
            H13 = H13 + DGDAi * DGDTi
            H22 = H22 + DGDFi * DGDFi
            H23 = H23 + DGDFi * DGDTi
            H33 = H33 + DGDTi * DGDTi

            F1 = F1 + Gi * DGDAi
            F2 = F2 + Gi * DGDFi
            F3 = F3 + Gi * DGDTi

        # Cholesly decomposition
        for adjust in range(1, AdjustMax):

            G11 = math.sqrt(H11 + VarLambda)
            G12 = H12 / G11
            G13 = H13 / G11
            G22 = math.sqrt(H22 + VarLambda - G12 * G12)
            G23 = (H23 - G12 * G13) / G22
            G33 = math.sqrt(H33 + VarLambda - G13 * G13 - G23 * G23)

            D1 = F1 / G11
            D2 = (F2 - G12 * D1) / G22
            D3 = (F3 - G13 * D1 - G23 * D2) / G33

            dT = D3 / G33
            dF = (D2 - G23 * dT) / G22
            dA = (D1 - G12 * dF - G13 * dT) / G11

            # Update the parameters
            Anew = Aold - dA
            Fnew = Fold - dF
            Tnew = Told - dT

            if 1 + 4 * Anew * Fnew < epsilon and VarLambda > 1:

                Xshift = Xshift + dX
                Yshift = Yshift + dY

                H = math.sqrt(1 + 4 * Aold * Fold)
                aTemp = -H * math.cos(Told) / (Aold + Aold) + dX
                bTemp = -H * math.sin(Told) / (Aold + Aold) + dY
                rTemp = 1 / abs(Aold + Aold)

                Anew = 1 / (rTemp + rTemp)
                aabb = aTemp * aTemp + bTemp * bTemp
                Fnew = (aabb - rTemp * rTemp) * Anew
                Tnew = math.acos(-aTemp / math.sqrt(aabb))

                if bTemp > 0:
                    Tnew = 2 * math.pi - Tnew

                VarNew = VarOld
                break

            if 1 + 4 * Anew * Fnew < epsilon:
                VarLambda = VarLambda * factorUp
                continue

            DD = 1 + 4 * Anew * Fnew
            D = math.sqrt(DD)
            CT = math.cos(Tnew)
            ST = math.sin(Tnew)

            # Initialize the sum
            GG = 0

            for i in range(0, n):

                Xi = X[i] + Xshift
                Yi = Y[i] + Yshift
                Zi = Xi * Xi + Yi * Yi
                Ui = Xi * CT + Yi * ST

                ADF = Anew * Zi + D * Ui + Fnew
                SQ = math.sqrt(4 * Anew * ADF + 1)
                DEN = SQ + 1
                Gi = 2 * ADF / DEN
                GG = GG + Gi * Gi

            VarNew = GG / (n - 3)

            H = math.sqrt(1 + 4 * Anew * Fnew)
            anew = -H * math.cos(Tnew) / (Anew + Anew) - Xshift
            bnew = -H * math.sin(Tnew) / (Anew + Anew) - Yshift
            Rnew = 1 / abs(Anew + Anew)

            if VarNew <= VarOld:
                progress = (abs(anew - aold) + abs(bnew - bold) + abs(Rnew - Rold)) / (Rnew + Rold)

                if progress < epsilon:
                    Aold = Anew
                    Fold = Fnew
                    Told = Tnew
                    VarOld = VarNew
                    finish = 1
                    break

                VarLambda = VarLambda * factorDown
                break

            else:  # no improvement
                VarLambda = VarLambda * factorUp
                continue

        if finish == 1:
            break

    H = math.sqrt(1 + 4 * Aold * Fold)
    result_a = -H * math.cos(Told) / (Aold + Aold) - Xshift
    result_b = -H * math.sin(Told) / (Aold + Aold) - Yshift
    result_r = 1 / abs(Aold + Aold)

    return(result_a, result_b, result_r)


def get_SSE(a, b, r, x, y):
    """
    Determine the sum of the squares of the errors (SSE) for a circle fit, helper function for calculating the curvature

    input: a, b, r, x, y.  circle center, radius, xpts, ypts
    output: SSE
    """

    SSE = 0

    for i in range(len(y)):
        xi = x[i]
        yi = y[i]
        v = (math.sqrt((xi - a)**2 + (yi - b)**2) - r)**2
        SSE += v

    return(SSE)
