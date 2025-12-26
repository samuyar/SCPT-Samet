from math import *


class Gleichgewicht:


    def __init__(self, As, Es, a, d, q, F, Fcc, Fct, FaiPa, FaiOr, Ffpz, Vuncr, Vda, beta1, beta2, z, zct, zai, zfpz,
                 zuncr, lambdaCS, x1, Ne, fy, fu, epsilonY, epsilonU, h3, scr, b, h, x0, n, ds, fct, fcm):
        """
        Inititierung des Gleichgewichts, Bestimmung/Berechnung der Einwirkungen und der Stahlzugkraft
        :param As
        :param Es
        :param a
        :param d
        :param q
        :param F
        :param Fcc
        :param Fct
        :param FaiPa
        :param FaiOr
        :param Ffpz
        :param Vuncr
        :param Vda
        :param beta1
        :param beta2
        :param z
        :param zct
        :param zai
        :param zfpz
        :param zuncr
        :param lambdaCS
        :param x1
        """

        self.Ve = - FaiOr * cos(beta2) + FaiPa * sin(beta2) + Ffpz * cos(
            beta1) + Vuncr + Vda  # Bestimmung der Einwirkenden Querkraft, keine Fallunterscheidung in Abhängigkeit von x1


        # Entscheidungsfunktion für die Bestimmung von Me (Einwirkendes Moment infolge des statischen Systems und der Belastung), keine Fallunterscheidung
        spanType = 1  # Belastungszustands derzeit festgelegt
        if spanType == 1:
            self.Me = self.Ve * lambdaCS * d  # type 1
        elif spanType == 2:
            self.Me = self.Ve * lambdaCS * d * (1 + 1 / 2 * (lambdaCS * d) / (a - lambdaCS * d))  # type 2
        elif spanType == 3:
            self.Me = self.Ve * lambdaCS * d + (a ** 2 + (lambdaCS * d) ** 2) / 2 * q  # type 3
        elif spanType == 4:
            self.Me = self.Ve * lambdaCS * d + (a ** 2 + (lambdaCS * d) ** 2) / 2 * q + a * F  # type 4

        self.Fs = self.Me / z + FaiPa * cos(beta2) + FaiOr * (sin(beta2) - zai / z) + Ffpz * (
                    zfpz / z - sin(beta1)) + Vuncr * zuncr / z + Ne * (1 - d / (2 * z)) + Fct * (
                              zct / z - 1)  # Stahlzugkraft

        # Vereinfachte Berechnung der Stahldehnung, später Berücksichtigung der verschiedener Einflüsse/Anteile auf/der Dehnung ()
        #self.epsilonS = self.Fs / (As * Es)

        # Vereinfachte Berechnung der Stahldehnung, bilinearer Ansatz
        self.sigmaS = self.Fs / As
        self.scrdelam = h3

        if self.sigmaS < fy:
            self.epsilonSdelam = self.Fs / (As * Es)

        else:
            self.epsilonSdelam = (fy / Es) + (self.sigmaS - fy) / ((fu - fy) / (epsilonU - epsilonY))

        self.epsilonS = self.epsilonSdelam

        # Berechnung der Stahldehnung unter Berücksichtigung der Delamination und des Tension Stiffening
        # Formeln nur gültig für große lambdacs (mit  d*lambdacs > scr)
        # für kleine lambdacs könnte der horizontale Riss das Auflager erreichen (siehe [Cla20] S.13)

        # Anteil aus Delamination
        #self.sigmaS = self.Fs / As
        #self.scrdelam = h3

        #if self.sigmaS < fy:
        #    self.epsilonSdelam = self.Fs / (As * Es)

        #else:
        #    self.epsilonSdelam = (fy / Es) + (self.sigmaS - fy) / ((fu - fy) / (epsilonU - epsilonY))

        # Anteil aus Tension Stiffening (mit linearer Annäherung an den Ansatz von Bentz)
        #self.hceff = min(2.5 * (h - d), (h - x0) / 3)
        #self.Aceff = self.hceff * b
        #self.M = self.Aceff / (n * ds * pi)

        #self.fyts = epsilonY * Es + fct / (1 + sqrt(3.6 * self.M * epsilonY))

        #if self.sigmaS < fy:
        #    self.epsilonSts = (self.sigmaS - fct) / ((self.fyts - fct) / epsilonY)
        #    if self.epsilonSts < 0:
        #        self.epsilonSts = 0  # ODER mit max(..., 0) 2 Zeilen zuvor umsetzen

        #else:
        #    self.epsilonSts = (fy / Es) + (self.sigmaS - fy) / ((fu - fy) / (epsilonU - epsilonY))

        # Berechnung der gesamten Stahldehnung
        #if self.scrdelam > scr:
        #    self.epsilonS = self.epsilonSdelam

        #else:
        #    self.epsilonS = (self.epsilonSts * (scr - self.scrdelam) + self.epsilonSdelam * self.scrdelam) / scr
