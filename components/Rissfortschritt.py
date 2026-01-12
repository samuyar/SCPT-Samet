from math import *


class Rissfortschritt:
    """
    # Berechnung des Rissfortschritts im Bauteil (x1) und der dazugörigen Spannungen
    # [mm, N/mm², Bogenmaß]
    """

    def __init__(self, d, scr, x0, y1, y2, beta, phi, fct, fcm, Ec, sigmaZ0, interlayer_status):
        """
        Initialisierung der Berechnung des Rissfortschritts im Bauteil (=x1) und der dazugörigen Spannungen
        :param d       statische Nutzhöhe - [mm]
        :param scr     Breite des Betonzahns - [mm]
        :param x0      Druckzonenhöhe - [mm]
        :param y1      vertikaler Längenanteil des Abschnitt Bs des Risses - [mm]
        :param y2      vertikaler Längenanteil des Abschnitt As des Risses - [mm]
        :param beta    Rissneigungswinkel an der Rissspitze - [Bogenmaß]
        :param phi     Stabverdrehungswinkel - [Bogenmaß]
        :param fct     Betonzugfestigkeit - [N/mm²]
        :param fcm     mittlere Betondruckfestigkeit - [N/mm²]
        :param Ec      mittlerer E-Modul - [N/mm²]
        :param sigmaZ0 Druckkraft aus Zahnbiegung - [N/mm²]
        """

        #Fügt den Boolean von interlayer_status hinzu
        self.interlayer_status = interlayer_status

        # wenn interlayer_status == True, dann Interlayer Werte einfügen (erstmal Reduzierung um 10%)
        if self.interlayer_status:
            fcm = 0.9 * fcm
            fct = 0.9 * fct
            Ec = 0.9 * Ec


        # Entscheidungfunktion zur Bestimmung von sigma1 in Abhängigkeit von der Lage der Rissspitze (Kupfer'sches Bruchkiterium)
        self.sigma1 = ((1 + 0.8 * sigmaZ0 / fcm * (1 + 1 / (tan(beta) ** 2))) / (
                    1 + 0.8 * fct / (fcm * (tan(beta) ** 2)))) * fct  # Hauptspannung 1 - [N/mm²]
        if self.sigma1 > fct:
            self.sigma1 = fct

        self.sigmaX0 = self.sigma1 * (1 - 1 / (tan(beta)) ** 2) + sigmaZ0 / (tan(
            beta) ** 2)  # Berechnung der Spannungen an der Rissspitze in X-Richtung mittels der Hauptspannung sigma1

        self.x1 = scr * self.sigmaX0 / (phi * Ec)

        self.tau0 = (self.sigma1 - sigmaZ0) / tan(beta)  # Schubspannung an Rissspitze - [N/mm²]

        self.sigma2 = -self.tau0 / tan(beta) + sigmaZ0  # Hauptspannung 2 - [N/mm²]

        self.alpha = (self.tau0 * tan(beta) + sigmaZ0) / (
                    - self.tau0 / tan(beta) + sigmaZ0)  # Verhältnis der Hauptspannungen zueinander

        # Kupfer'sches Bruchkriterium für weitere Berechnung - [N/mm²]

    #        self.sigmaBiFrac = None;
    #        if self.alpha >= 0:
    #            self.sigmaBiFrac = fct;
    #        else:
    #            self.sigmaBiFrac = (1 + 0.8 * self.sigma2/fcm) * fct;

    # Einbau von Gettern für den späteren Gebrauch

    def getX1(self):
        return self.x1

    def getSigmaX0(self):
        return self.sigmaX0

    def getSigma1(self):
        return self.sigma1

    def getSigma2(self):
        return self.sigma2

    def getTau0(self):
        return self.tau0

    def getAlpha(self):
        return self.alpha

    def get_interlayer_status(self):
        return self.interlayer_status
