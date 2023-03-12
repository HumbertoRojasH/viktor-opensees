import pandas as pd
import numpy as np
from openseespy.opensees import *
from numpy import zeros
from viktor import Color
from viktor.geometry import Sphere, Point, Material, Group, Line, Polyline
from viktor.utils import memoize
import plotly.graph_objects as go


@memoize
def GeoModel(dx, dy, h, nx, ny, nz):
    from numpy import zeros
    Lx, Ly, Lz = dx * nx, dy * ny, h * nz
    NN = (nx + 1) * (ny + 1) * (nz + 1)
    Nodes = zeros((NN, 5))
    # Creando los nodos y asignando coordenadas
    # Creating the nodes and assigning coordinates
    c = 0
    for i in range(nz + 1):
        for j in range(ny + 1):
            for k in range(nx + 1):
                if k == nx and j != ny and j != 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.50]
                elif k != nx and j == ny and k != 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.50]
                elif k == 0 and j != ny and j != 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.50]
                elif k != nx and j == 0 and k != 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.50]
                elif k == nx and j == ny:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.25]
                elif k == 0 and j == 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.25]
                elif k == nx and j == 0:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.25]
                elif k == 0 and j == ny:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 0.25]
                else:
                    Nodes[c] = [c, k * dx, j * dy, i * h, 1.00]
                c = c + 1
    Nodes[:(nx + 1) * (ny + 1), 4] = 0

    NE = (nx * (ny + 1) + ny * (nx + 1) + (nx + 1) * (ny + 1)) * nz
    Elems = zeros((NE, 5), 'int')
    # Creando las conexiones de los elementos verticales
    # Creating the connections of the vertical elements
    c = 0
    for i in range(nz):
        for j in range(ny + 1):
            for k in range(nx + 1):
                Elems[c] = [c, c, c + (nx + 1) * (ny + 1), 1, 0]
                c = c + 1
    # Creando las conexiones de los elementos horizontales
    # Creating the connections of the horizontal elements
    m = (nx + 1) * (ny + 1)
    for i in range(nz):
        for j in range(ny + 1):
            for k in range(nx):
                Elems[c] = [c, m, m + 1, 2, 0]
                m = m + 1
                c = c + 1
            m = m + 1
    # Creando las conexiones de los elementos horizontales
    # Creating the connections of the horizontal elements
    n = 0
    for i in range(nz):
        n = n + (nx + 1) * (ny + 1)
        for j in range(nx + 1):
            for k in range(ny):
                Elems[c] = [c, j + k * (nx + 1) + n, j + nx + 1 + k * (nx + 1) + n, 2, 1]
                c = c + 1
    # Creando centro de diafragmas
    # Creating center of diaphragms
    Diap = zeros((nz, 4))
    for i in range(nz):
        Diap[i] = [i + 1000, Lx / 2.0, Ly / 2.0, h * (i + 1)]

    return Nodes, Elems, Diap


def espectro_E030(T, Z=0.45, U=1.5, S=1.0, Tp=0.4, Tl=2.5, R=1):
    # Definimos el factor de cortante basal de acuerdo a la norma peruana E.030
    # We define the basal shear factor according to the Peruvian standard E.030
    from numpy import zeros
    n = len(T)
    E030 = zeros(n)
    for i in range(n):
        if T[i] >= 0 and T[i] < 0.2 * Tp:
            E030[i] = 2.5  # 1+7.5*T[i]/Tp
        elif T[i] >= 0.2 * Tp and T[i] < Tp:
            E030[i] = 2.5
        elif T[i] >= Tp and T[i] < Tl:
            E030[i] = 2.5 * (Tp / T[i])
        elif T[i] >= Tl:
            E030[i] = 2.5 * (Tp * Tl / T[i] ** 2)
        else:
            print("El periodo no puede ser negativo!")
    return E030 * Z * U * S / R


def get_static_loads(coef, p, h, T):
    # Obteniendo las fuerzas laterales en cada entrepiso de acuerdo a la norma peruana E.030
    # Obtaining the lateral forces in each mezzanine according to the Peruvian standard E.030
    from numpy import zeros
    n = len(h)
    V = coef * sum(p)
    F = zeros(n)
    #
    if T > 0.0 and T <= 0.5:
        k = 1.0
    elif T > 0.5:
        k = 0.75 + 0.5 * T
    else:
        print('El periodo es negativo!')
    #
    div = 0.
    for i in range(n):
        div = div + p[i] * h[i] ** k
    #
    for i in range(n):
        F[i] = p[i] * h[i] ** k / div * V
    return F, k


@memoize
def analysis_static(Nodes, Elems, Diap, params):
    s1 = params.s1
    s2 = params.s2
    # Unidades
    # units

    # Unidades Base
    # Base Units
    m = 1
    kg = 1
    s = 1
    # Otras Unidades
    # Other Units
    cm = 0.01 * m
    N = kg * m / s ** 2
    kN = 1000 * N
    kgf = 9.80665 * N
    tonf = 1000 * kgf
    Pa = N / m ** 2
    # Constantes Físicas
    # Physical Constants
    g = 9.80665 * m / s ** 2

    # Propiedades del Concreto
    # Concrete Properties
    fc = s1.fc * kg / cm ** 2
    E = 151 * fc ** 0.5 * kgf / cm ** 2
    G = 0.5 * E / (1 + 0.2)

    # Propiedades de la Viga 
    # Beam Properties
    Av = s1.vb * s1.vh
    Izv = s1.vb * s1.vh ** 3 / 12
    Iyv = s1.vb ** 3 * s1.vh / 12
    aa, bb = max(s1.vb, s1.vh), min(s1.vb, s1.vh)
    β = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa) ** 4 / 12)
    Jxxv = β * bb ** 3 * aa

    # Propiedades de la Columna
    # Column Properties
    Ac = s1.cb * s1.ch
    Izc = s1.cb * s1.ch ** 3 / 12
    Iyc = s1.cb ** 3 * s1.ch / 12
    aa, bb = max(s1.cb, s1.ch), min(s1.cb, s1.ch)
    β = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa) ** 4 / 12)
    Jxxc = β * bb ** 3 * aa

    wipe()
    model('basic', '-ndm', 3, '-ndf', 6)
    RigidDiaphragm = 'ON'

    # Creamos los nodos
    # Creating nodes
    for Ni in Nodes:
        node(int(Ni[0]), *Ni[1:4])

    # Definimos diafragmas rígidos
    # Defining rigid diaphragms
    if RigidDiaphragm == 'ON':
        # perpendicular al plano del diafragma
        # perpendicular to the plane of the diaphragm
        dirDia = 3
        for Nd in Diap:
            node(int(Nd[0]), *Nd[1:4])
            fix(int(Nd[0]), *[0, 0, 1, 1, 1, 0])
            NodesDi = []
            for Ni in Nodes:
                if Ni[3] == Nd[3]:
                    NodesDi.append(int(Ni[0]))
            rigidDiaphragm(dirDia, int(Nd[0]), *NodesDi)

    # Restricciones
    # Restraints
    fixZ(0.0, *[1, 1, 1, 1, 1, 1], '-tol', 1e-6)

    geomTransf('PDelta', 1, *[1, 0, 0])
    geomTransf('Linear', 2, *[1, -1, 0])

    # Creamos los elementos
    # Creating elements
    for Ele in Elems:
        if int(Ele[3]) == 1:  # 1 Columna # 1 Columns
            element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Ac, E, G, Jxxc, Iyc, Izc, int(Ele[3]),
                    '-mass', s1.d * Ac * 10 ** -8)
        else:  # 2 Viga # 2 Beam
            element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Av, E, G, Jxxv, Iyv, Izv, int(Ele[3]),
                    '-mass', s1.d * Av * 10 ** -8)

    wLive = s2.sec2.cv
    wLosa = s2.sec2.cl
    wAcab = s2.sec2.ca
    wTabi = s2.sec2.ct
    wTotal = 1.0 * (wLosa + wAcab + wTabi) + 0.25 * wLive
    Carga = wTotal * s1.ex * s1.ey * m ** 2

    for Ni in Nodes:
        mass(int(Ni[0]), Ni[4] * Carga, Ni[4] * Carga, 0.0)

    Nmodes = 3 * s1.nz
    vals = eigen(Nmodes)
    Tmodes = np.zeros(len(vals))
    for i in range(Nmodes):
        Tmodes[i] = 2 * np.pi / vals[i] ** 0.5

    # Realizamos un análisis para obtener la matriz de Masas
    # We perform an analysis to obtain the Mass Matrix
    wipeAnalysis()
    system('FullGeneral')
    numberer("Plain")
    constraints('Transformation')
    algorithm('Linear')
    analysis('Transient')
    integrator('GimmeMCK', 1.0, 0.0, 0.0)
    analyze(1, 0.0)

    # Obtenemos la matriz de Masas
    # We get the Mass Matrix
    N = systemSize()  # Número de Grados de Libertad # Number of Degrees of Freedom
    Mmatrix = printA('-ret')
    Mmatrix = np.array(Mmatrix).reshape((N, N))
    MF = Mmatrix[-3 * s1.nz:, -3 * s1.nz:]
    H = np.arange(1, s1.nz + 1) * s1.ez
    P = sum(MF[0::3, 0::3]) * g

    E030 = espectro_E030(Tmodes, s2.sec3.z, s2.sec3.u, s2.sec3.s, s2.sec3.tp, s2.sec3.tl, s2.sec3.r)
    F, k = get_static_loads(E030[0], P, H, Tmodes[0])

    if s2.sec1.dir == "Sismo X-X":
        timeSeries('Linear', 1)
        pattern('Plain', 1, 1)
        Le = s1.ny * s1.ey * 0.05
        for i in range(s1.nz):
            load(int(Diap[i][0]), F[i], 0., 0., 0., 0., F[i] * Le)

    elif s2.sec1.dir == "Sismo Y-Y":
        timeSeries('Linear', 1)
        pattern('Plain', 1, 1)
        Le = s1.nx * s1.ex * 0.05
        for i in range(s1.nz):
            load(int(Diap[i][0]), 0., F[i], 0., 0., 0., F[i] * Le)

    wipeAnalysis()
    constraints('Transformation')
    numberer('Plain')
    system('FullGeneral')
    algorithm('Linear')
    integrator('LoadControl', 1)
    analysis('Static')
    analyze(1)

    VS = np.cumsum(F[::-1])[::-1]
    df3 = pd.DataFrame(
        columns=['Nivel', 'H(m)', 'V(kN)', 'Ux(cm)', 'Uy(cm)', '0.75R*Ux(cm)', '0.75R*Uy(cm)', 'DriftX(‰)',
                 'DriftY(‰)'])
    tempX, tempY = 0., 0.

    u = zeros(3 * len(Nodes))
    for i in range(len(Nodes)):
        u[3 * i] = nodeDisp(int(Nodes[i][0]), 1)
        u[3 * i + 1] = nodeDisp(int(Nodes[i][0]), 2)
        u[3 * i + 2] = nodeDisp(int(Nodes[i][0]), 3)

    for i in range(s1.nz):
        desX = nodeDisp(int(Diap[i][0]), 1)
        desY = nodeDisp(int(Diap[i][0]), 2)

        rotZ = nodeDisp(int(Diap[i][0]), 6)
        desX = desX + abs(rotZ * s1.ny * s1.ey / 2)
        desY = desY + abs(rotZ * s1.nx * s1.ex / 2)

        desXR, desYR = desX * 0.75 * s2.sec3.r, desY * 0.75 * s2.sec3.r
        driftX = 1000. * (desXR - tempX) / s1.ez
        driftY = 1000. * (desYR - tempY) / s1.ez
        tempX, tempY = desXR, desYR

        df3 = df3.append(
            {'Nivel': i + 1, 'H(m)': H[i], 'V(kN)': VS[i] / 1000, 'Ux(cm)': desX * 100, 'Uy(cm)': desY * 100,
             '0.75R*Ux(cm)': desXR * 100, '0.75R*Uy(cm)': desYR * 100,
             'DriftX(‰)': driftX, 'DriftY(‰)': driftY}, ignore_index=True)

    return df3.round(4), u


def getCombo(E030, MF, modo, Tmodes, NT, ni):
    # Definimos valores iniciales
    # We define initial values
    D_ABSx, D_RCSCx = np.zeros(NT), np.zeros(NT)
    Δ_ABSx, Δ_RCSCx = np.zeros(NT), np.zeros(NT)
    V_ABSx, V_RCSCx = np.zeros(NT), np.zeros(NT)
    D_ABSy, D_RCSCy = np.zeros(NT), np.zeros(NT)
    Δ_ABSy, Δ_RCSCy = np.zeros(NT), np.zeros(NT)
    V_ABSy, V_RCSCy = np.zeros(NT), np.zeros(NT)
    
    Ux, Uy, Rz = np.zeros(NT), np.zeros(NT), np.zeros(NT)
    Ux[0::3] = 1
    Uy[1::3] = 1
    Rz[2::3] = 1

    # Se realiza la Superpocisión Modal Espectral
    # Spectral Modal Superposition is performed
    for j in range(1, ni + 1):  # ni+1
        FPx = modo[j - 1].T @ MF @ Ux
        FPy = modo[j - 1].T @ MF @ Uy
        FPr = modo[j - 1].T @ MF @ Rz
        #
        Sa = E030[j - 1] * 9.80665
        Sd = Sa / (2 * np.pi / Tmodes[j - 1]) ** 2
        #
        respDX = Sd * FPx * modo[j - 1]
        respAX = Sa * FPx * MF @ modo[j - 1]
        D_ABSx = D_ABSx + abs(respDX)
        D_RCSCx = D_RCSCx + (respDX) ** 2
        respDX[3:] = respDX[3:] - respDX[:-3]
        Δ_ABSx = Δ_ABSx + abs(respDX)
        Δ_RCSCx = Δ_RCSCx + (respDX) ** 2
        V_ABSx = V_ABSx + abs(np.cumsum(respAX[::-1])[::-1])
        V_RCSCx = V_RCSCx + (np.cumsum(respAX[::-1])[::-1]) ** 2
        #
        respDY = Sd * FPy * modo[j - 1]
        respAY = Sa * FPy * MF @ modo[j - 1]
        D_ABSy = D_ABSy + abs(respDY)
        D_RCSCy = D_RCSCy + (respDY) ** 2
        respDY[3:] = respDY[3:] - respDY[:-3]
        Δ_ABSy = Δ_ABSy + abs(respDY)
        Δ_RCSCy = Δ_RCSCy + (respDY) ** 2
        V_ABSy = V_ABSy + abs(np.cumsum(respAY[::-1])[::-1])
        V_RCSCy = V_RCSCy + (np.cumsum(respAY[::-1])[::-1]) ** 2

    # Se realiza la combinación 25%ABS + 75%RCSC
    # 25%ABS + 75%RCSC combination is performed
    D_RCSCx = D_RCSCx ** 0.5
    Δ_RCSCx = Δ_RCSCx ** 0.5
    V_RCSCx = V_RCSCx ** 0.5
    DDx = 0.25 * D_ABSx + 0.75 * D_RCSCx
    ΔDx = 0.25 * Δ_ABSx + 0.75 * Δ_RCSCx
    VDx = 0.25 * V_ABSx + 0.75 * V_RCSCx
    #
    D_RCSCy = D_RCSCy ** 0.5
    Δ_RCSCy = Δ_RCSCy ** 0.5
    V_RCSCy = V_RCSCy ** 0.5
    DDy = 0.25 * D_ABSy + 0.75 * D_RCSCy
    ΔDy = 0.25 * Δ_ABSy + 0.75 * Δ_RCSCy
    VDy = 0.25 * V_ABSy + 0.75 * V_RCSCy

    df = pd.DataFrame(columns=['Nivel', 'VDx(kN)', 'VDy(kN)', 'UDx(cm)', 'UDy(cm)'])
    for i in range(int(NT / 3)):
        df = df.append({'Nivel': i + 1, 'VDx(kN)': VDx[0::3][i] / 1000,
                        'VDy(kN)': VDy[1::3][i] / 1000, 'UDx(cm)': DDx[0::3][i] * 100,
                        'UDy(cm)': DDy[1::3][i] * 100}, ignore_index=True)

    return DDx, ΔDx, VDx, DDy, ΔDy, VDy, df


@memoize
def analysis_modal(Nodes, Elems, Diap, params):
    s1 = params.s1
    s2 = params.s2
    # Unidades
    # Units

    # Unidades Base
    # Base Units
    m = 1
    kg = 1
    s = 1
    # Otras Unidades
    # Other Units
    cm = 0.01 * m
    N = kg * m / s ** 2
    kN = 1000 * N
    kgf = 9.80665 * N
    tonf = 1000 * kgf
    Pa = N / m ** 2
    # Constantes Físicas
    # Physical Constants
    g = 9.80665 * m / s ** 2

    # Propiedades del Concreto
    # Concrete Properties
    fc = s1.fc*kg/cm**2
    E = 151*fc**0.5*kgf/cm**2
    G = 0.5 * E / (1 + 0.2)

    # Propiedades de la Viga
    # Beam Properties
    Av = s1.vb * s1.vh
    Izv = s1.vb * s1.vh ** 3 / 12
    Iyv = s1.vb ** 3 * s1.vh / 12
    aa, bb = max(s1.vb, s1.vh), min(s1.vb, s1.vh)
    β = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa) ** 4 / 12)
    Jxxv = β * bb ** 3 * aa

    # Propiedades de la Columna
    # Column Properties
    Ac = s1.cb * s1.ch
    Izc = s1.cb * s1.ch ** 3 / 12
    Iyc = s1.cb ** 3 * s1.ch / 12
    aa, bb = max(s1.cb, s1.ch), min(s1.cb, s1.ch)
    β = 1 / 3 - 0.21 * bb / aa * (1 - (bb / aa) ** 4 / 12)
    Jxxc = β * bb ** 3 * aa

    wipe()
    model('basic', '-ndm', 3, '-ndf', 6)
    RigidDiaphragm = 'ON'

    # Creamos los nodos
    # Create nodes
    for Ni in Nodes:
        node(int(Ni[0]), *Ni[1:4])

    # Definimos diafragmas rígidos
    # Define rigid diaphragms
    if RigidDiaphragm == 'ON':
        dirDia = 3  # perpendicular al plano del diafragma # perpendicular to the plane of the diaphragm
        for Nd in Diap:
            node(int(Nd[0]), *Nd[1:4])
            fix(int(Nd[0]), *[0, 0, 1, 1, 1, 0])
            NodesDi = []
            for Ni in Nodes:
                if Ni[3] == Nd[3]:
                    NodesDi.append(int(Ni[0]))
            rigidDiaphragm(dirDia, int(Nd[0]), *NodesDi)

    # Restricciones
    # Restraints
    fixZ(0.0, *[1, 1, 1, 1, 1, 1], '-tol', 1e-6)

    geomTransf('PDelta', 1, *[1, 0, 0])
    geomTransf('Linear', 2, *[1, -1, 0])

    # Creamos los elementos
    # Create elements
    for Ele in Elems:
        if int(Ele[3]) == 1:  # 1 Columna # 1 Columns
            element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Ac, E, G, Jxxc, Iyc, Izc, int(Ele[3]),
                    '-mass', s1.d * Ac)
        else:  # 2 Viga # 2 Beam
            element('elasticBeamColumn', int(Ele[0]), int(Ele[1]), int(Ele[2]), Av, E, G, Jxxv, Iyv, Izv, int(Ele[3]),
                    '-mass', s1.d * Av)

    wLive = s2.sec2.cv
    wLosa = s2.sec2.cl
    wAcab = s2.sec2.ca
    wTabi = s2.sec2.ct
    wTotal = 1.0 * (wLosa + wAcab + wTabi) + 0.25 * wLive
    Carga = wTotal * s1.ex * s1.ey

    for Ni in Nodes:
        mass(int(Ni[0]), Ni[4] * Carga, Ni[4] * Carga, 0.0)

    nz = s1.nz
    Nmodes = 3*nz
    vals = eigen(Nmodes)
    Tmodes = np.zeros(len(vals))
    for i in range(Nmodes):
        Tmodes[i] = 2*np.pi/vals[i]**0.5

    # Realizamos un análisis para obtener la matriz de Masas
    # We perform an analysis to obtain the Mass Matrix
    wipeAnalysis()
    system('FullGeneral')
    numberer("Plain")
    constraints('Transformation')
    algorithm('Linear')
    analysis('Transient')
    integrator('GimmeMCK', 1.0, 0.0, 0.0)
    analyze(1, 0.0)

    # Obtenemos la matriz de Masas
    # We get the Mass Matrix
    N = systemSize()  # Número de Grados de Libertad # Number of Degrees of Freedom
    Mmatrix = printA('-ret')
    Mmatrix = np.array(Mmatrix).reshape((N, N))
    MF = Mmatrix[-3 * nz:, -3 * nz:]
    H = np.arange(1, nz + 1) * s1.ez
    P = sum(MF[0::3, 0::3]) * 9.80665  # Peso por nivel # weight per level
    E030 = espectro_E030(Tmodes, s2.sec3.z, s2.sec3.u, s2.sec3.s, s2.sec3.tp, s2.sec3.tl, s2.sec3.r)
    F, k = get_static_loads(E030[0], P, H, Tmodes[0])
    VS = np.cumsum(F[::-1])[::-1]

    Tags = getNodeTags()
    modo = np.zeros((Nmodes, 3 * nz))
    for j in range(1, Nmodes + 1):
        ind = 0
        for i in Tags[-nz:]:
            temp = nodeEigenvector(i, j)
            modo[j - 1, [ind, ind + 1, ind + 2]] = temp[0], temp[1], temp[-1]
            ind = ind + 3

    # Definimos valores iniciales
    # Define initial values
    Ux, Uy, Rz = np.zeros(3 * nz), np.zeros(3 * nz), np.zeros(3 * nz)
    Ux[0::3] = 1
    Uy[1::3] = 1
    Rz[2::3] = 1
    SUMx, SUMy, SUMr = 0., 0., 0.
    ni = 0

    Mx = sum(sum(MF[0::3, 0::3]))
    My = sum(sum(MF[1::3, 1::3]))
    Mr = sum(sum(MF[2::3, 2::3]))

    df3 = pd.DataFrame(columns=['Modo', 'T(s)', 'SumUx', 'SumUy', 'SumRz'])
    for j in range(1, Nmodes + 1):
        FPx = modo[j - 1].T @ MF @ Ux
        FPy = modo[j - 1].T @ MF @ Uy
        FPr = modo[j - 1].T @ MF @ Rz
        FPRx = FPx ** 2 / Mx
        FPRy = FPy ** 2 / My
        FPRr = FPr ** 2 / Mr
        SUMx = SUMx + FPRx
        SUMy = SUMy + FPRy
        SUMr = SUMr + FPRr
        #
        if min(SUMx, SUMy, SUMr) >= 0.90 and ni == 0:
            ni = j
        df3 = df3.append({'Modo': j, 'T(s)': Tmodes[j - 1], 'SumUx': SUMx,
                          'SumUy': SUMy, 'SumRz': SUMr}, ignore_index=True)

    DDx, ΔDx, VDx, DDy, ΔDy, VDy, df4 = getCombo(E030, MF, modo, Tmodes, 3 * nz, ni)
    df4 = df4.astype({'Nivel': int})

    # Escalamiento de los resultados del análisis dinámico
    # Scaling of the results of the dynamic analysis
    if VDx[0::3][0] < 0.80 * VS[0]:
        FSx = 0.80 * VS[0] / VDx[0::3][0]
    else:
        FSx = 1.
        msjx = 'NO es necesario escalar en X'

    if VDy[1::3][0] < 0.80 * VS[0]:
        FSy = 0.80 * VS[0] / VDy[1::3][0]
    else:
        FSy = 1.

    H = np.arange(1, nz + 1) * s1.ez
    df5 = pd.DataFrame(columns=['Nivel', 'H(m)', 'Vx(kN)', 'Vy(kN)', 'Ux(cm)', 'Uy(cm)', 'Δx(‰)', 'Δy(‰)'])
    for i in range(nz):
        Δx = 0.75 * s2.sec3.r * ΔDx[0::3][i] / s1.ez
        Δy = 0.75 * s2.sec3.r * ΔDy[1::3][i] / s1.ez
        #
        df5 = df5.append({'Nivel': i + 1, 'H(m)': H[i], 'Vx(kN)': FSx * VDx[0::3][i] / 1000,
                          'Vy(kN)': FSy * VDy[1::3][i] / 1000, 'Ux(cm)': 0.75 * s2.sec3.r * DDx[0::3][i] * 100,
                          'Uy(cm)': 0.75 * s2.sec3.r * DDy[1::3][i] * 100, 'Δx(‰)': Δx * 1000, 'Δy(‰)': Δy * 1000},
                         ignore_index=True)
    df5 = df5.astype({'Nivel': int})

    return df3.round(4), df5.round(4)


def displacement(nodes, U, scale):
    displacement = U.reshape(len(nodes), 3)
    deform = displacement + nodes
    deform_scale = displacement * scale + nodes
    return deform_scale


def plot_structure(nodes, lines, deform_scale=None):
    fig = go.Figure()

    if deform_scale is None:
        for i, (x, y) in enumerate(lines):
            fig.add_trace(go.Scatter3d(x=nodes[[x, y], 0], y=nodes[[x, y], 1], z=nodes[[x, y], 2], mode='lines+markers',
                                       line=dict(color='black', width=3), marker=dict(size=3, color='black')))

    else:
        for i, (x, y) in enumerate(lines):
            fig.add_trace(go.Scatter3d(x=nodes[[x, y], 0], y=nodes[[x, y], 1], z=nodes[[x, y], 2], mode='lines',
                                       line=dict(color='gray', width=1, )))
        for i, (x, y) in enumerate(lines):
            fig.add_trace(go.Scatter3d(x=deform_scale[[x, y], 0], y=deform_scale[[x, y], 1], z=deform_scale[[x, y], 2],
                                       mode="markers + lines",
                                       line=dict(color='black', width=3), marker=dict(size=3, color='black')))
    fig.update_layout(
        autosize=False,
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0),
        scene=dict(
            xaxis=dict(showgrid=False, showbackground=False, zerolinecolor="white", showticklabels=False),
            yaxis=dict(showgrid=False, showbackground=False, zerolinecolor="white", showticklabels=False),
            zaxis=dict(showgrid=False, showbackground=False, zerolinecolor="white", showticklabels=False))
    )
    return fig


def plot_structure_geometry_view(nodes, lines, deform_scale=None):
    _lines = []
    _deformed_lines = []
    _nodes = []

    if deform_scale is None:
        for x, y in lines:
            _line = []
            for node in zip(nodes[[x, y], 0], nodes[[x, y], 1], nodes[[x, y], 2]):
                pnt = Point(*node)
                _nodes.append(Sphere(centre_point=pnt, radius=0.1, material=Material(name='black', color=Color.black())))
                _line.append(pnt)
            _lines.append(Polyline(_line))
    else:
        for x, y in lines:
            _line = []
            for node in zip(deform_scale[[x, y], 0], deform_scale[[x, y], 1], deform_scale[[x, y], 2]):
                pnt = Point(*node)
                _nodes.append(Sphere(centre_point=pnt, radius=0.1, material=Material(name='black', color=Color.black())))
                _line.append(pnt)
            _deformed_lines.append(Polyline(_line))
    return Group(_nodes + _lines + _deformed_lines)