from math import *
import scipy.optimize
from components.DefVar import *

class ParameterNeu:

    def __init__(self, lambdaCS, scr, d, b, h, phi, x0, x1, y1, y2, deltaY, beta, beta1, beta2, epsilonS, epsilonSdelam,
                 sigma1, Ne, Fs, Fcc, Ffpz, FaiOr, FaiPa, Vda, z, zai, zfpz, zct, zuncr, yvor, scrdelam, ds, ns, fcm):
        """
        Initiierung der Berechnung der neuen Paramter
        :param lambdaCS
        :param scr
        :param d
        :param b
        :param phi
        :param x0
        :param x1
        :param y1
        :param y2
        :param deltaY
        :param beta
        :param beta1
        :param beta2
        :param epsilonS
        :param sigmaBiFrac
        :param Ne
        :param Fs
        :param Fcc
        :param Ffpz
        :param FaiOr
        :param FaiPa
        :param Vda
        :param z
        :param zai
        :param zfpz
        :param zct
        """

        # Bestimmung der Rotation der Rissspitze um CR in Abhängigkeit von beta, Abbruchkriterium phi = konstant
        self.phiCalc = epsilonS * scr / (d - x0)

        # Korrektheit prüfen!
        # Bestimmung der vertikalen Spannung an der Rissspitze sigmaZ0 (Cantilever Action)
        # Bestimmen des Spannmomentes Mca
        self.n = lambdaCS * d / scr #eq 69
        #self.deltaFs = Fs * (self.n + 1) / self.n - Fs
        self.limitBS = (scr - scrdelam) * pi * (ds * ns) * 1.26 * sqrt(fcm / 20)  # eq 68, Verbundbedingung, MC2010, Berücksichtigung der Delamination
        self.deltaFs = min(Fs * (self.n + 1) / self.n - Fs, self.limitBS) #eq 68

        # Contribution epsTs from regions with intact bond controlled by tension stiffening effect of surrounding concrete
        # Proposal by Bentz (eq. (64)):
        hceff = min(2, 5 * (h - d), (h - x0) / 3)
        Aceff = hceff * b
        M = Aceff / (pi * ds * ns)

        def tsZero(eps: float) -> float: #Here eps = epsTs
           return eps * DefVar.Es + DefVar.fct / (1 + sqrt(3.6 * M * abs(eps))) - Fs / DefVar.As #eq. (64) (Fs/As must be equal to the right hand side of eq. (64))

        solution = scipy.optimize.newton(tsZero,0)    #Use Newton-Raphson iteration to solve for epsTs
        epsTs = solution          #First entry in solution array is the estimated location where the function is zero

        #Calculation of averaged value of reinforcement strain epsilonS (eq. (65)):

        if scrdelam <= scr:
            epsilonS = (epsTs * (scr - scrdelam) + epsilonSdelam * scrdelam) / scr #eq. (65)
        else:
            epsilonS = epsilonSdelam         #If delamination crack is present all along scr, then epsilonS = epsdelam (no bond)

        Mca = - Ffpz * scr * cos(beta1) + FaiOr * scr * cos(beta2) - FaiPa * scr * sin(beta2) - Vda * scr + self.deltaFs * (y2 + y1)
        Wc = b * scr ** 2 / 6  # Annahme für Widerstandsmoment
        self.sigmaZ0 = Mca / Wc  # Cantilever Action, Vorzeichenkonform?

        # Bestimmen des zusätzlichen Risswachstums deltaY
        self.deltaYneu = d - x0 - yvor - x1

        # Bestimmen des Neigungswinkels beta an der Rissspitze,  bisher nur für Span-Type 1 (Formel nicht allgemeingültig!)
        u = 2 / 3 * b * ((x0 + x1) + 1 / 2 * (x1 ** 2 + x0 * x1) / (x0 - x1)) * (sigma1 - self.sigmaZ0) * (zuncr + lambdaCS * d)
        v = Fcc * z + FaiOr * (zai + lambdaCS * d * cos(beta2)) - FaiPa * lambdaCS * d * sin(beta2) - Ffpz * (
                zfpz + lambdaCS * d * cos(beta1)) + Ne * d / 2 - Vda * lambdaCS * d - 1 / 2 * b * x1 * sigma1 * zct
        w = 1 / 2 * b * x1 * zct * (sigma1 - self.sigmaZ0)

        try:
            self.betaNeu = [atan(u / (2 * v) + sqrt((u / (2 * v)) ** 2 - w / v)),
                            atan(u / (2 * v) - sqrt((u / (2 * v)) ** 2 - w / v))]
            self.betaNeu.sort(
                reverse=True)  # sortieren der Elemente der Größe nach, Reihenfolge der Sortierung prüfen! Größter Wert zu erst
            self.betaNeu = self.betaNeu[0]
        except ValueError as error:
            self.betaNeu = 1.5
            #print(error)

