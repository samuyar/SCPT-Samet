from math import *


class Kinematik:
    """
    Berechnung der Rissöffnungen und -Verschiebungen
    [mm, N/mm², Bogenmaß]
    """

    def __init__(self, phi, scr, y1, y2, x0, x1, beta1, beta2):
        """
        :param phi     Stabverdrehungswinkel - [Bogenmaß]
        :param scr     Breite des Betonzahns nach Classen - [mm]
        :param y1      Vertikale Risslänge Abschnitt B - [mm]
        :param y2      Vertikale Risslänge Abschnitt A - [mm]
        :param x0      Höhe der Betondruckzone - [mm]
        :param x1      Abstand Center of Rotation zu Rissspitze - [mm]
        :param beta1   Rissneigung Abschnitt B - [Bogenmaß]
        :param beta2   Rissneigung Abschnitt A - [Bogenmaß]
        """
        self.l2 = y2/sin(beta2)                # Risslänge [mm]

        self.wfpz = phi * (y1 / sin(beta1) + sin(beta1) * x1)   # Rissöffnung unten Abschnitt B, ! Einschränkung für x1

        self.r1 = abs(y1/sin(beta1) * sin(beta2 - beta1) - x1 * cos(beta2))  # Radius Abschnitt B, Absoluter Wert zu benutzen?
        self.r2 = abs(y1/sin(beta1) * cos(beta2 - beta1) + x1 * sin(beta2))  # Radius Abschnitt A, absoluter Wert zu benutzen?

        self.delta = self.r1 * phi                             # Rissgleiten? Über den unteren Riss konstant
        self.deltak = phi * (y1/tan(beta1) + y2/tan(beta2))    # Vertikale Rissuferverschiebung

        self.epsilonCr = phi/scr * x1          # Horizontale Dehnung der die Rissspitze berührenden Faser

        # Prüfe: Werden nachfolgende Dehnungen in folgenden Berechnungen benötigt?
        self.epsilonTop = - phi/scr * x0       # Horizontale Dehnung der oberen Faser der Druckzone

        # Noch unklar, wird in späterem Programm verwendet (für Bestimmung der Kräfte in Konstitutiv (5b)):
        # self.wtop_0 = self.wfpz * 0
        # self.wtop_1 = self.wfpz * 1
        # wbot siehe Konstitutiv.py

    # Definition von Gettern für den Vereinfachten Gebrauch
    def getWfpz(self):
        return self.wfpz

    def getR1(self):
        return self.r1

    def getR2(self):
        return self.r2

    def getDelta(self):
        return self.delta

    def getDeltak(self):
        return self.deltak

    def getEpsilonCr(self):
        return self.epsilonCr

    def getEpsilonTop(self):
        return self.epsilonTop
