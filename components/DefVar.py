from math import *


class DefVar:
    """
    Alle Angaben in [mm],[N/mm²]
    Bisherige Vereinfachung: Einlagige Bewehrung
    """
    # Angaben zum System, Geometrische Grunddaten, Materialkennwerte, Angaben zur Bewehrung und Verhältniswert des Risses
    lSpann = 1440  # Angabe der Spannweite (Mitte Auflager bis Mitte Auflager) [mm]

    #a = lSpann / 2  # Halbe Spannweite, evtl. später in Controller direkt oder in Modul Gleichgewicht berechnen?
    a = 720

    F = 0  # Einwirkende Einzellast

    q = 0  # Einwirkende Linienlast

    # Hier später evtl. automatische Auswahl des Shear-Span-Typs (evtl An)

    Ne = 0

    def getNe(self):
        return self.Ne

    d = 180  # statische Nutznöhe [mm]

    def getD(self):
        return self.d

    lambdaCS = 2.73  # Ort des erzeugten Risses

    b = 100  # Breite [mm]

    def getB(self):
        return self.b

    h = 210  # Höhe [mm]

    def getH(self):
        return self.h

    fcm = 54  # mittlere Betondruckfestigkeit [N/mm²]

    def getFcm(self):
        return self.fcm

    #fct = 0.3 * (fcm - 4) ** (
    #        2 / 3)  # Betonzugfestigkeit aus Betondruckfestigkeit berechnet [N/mm²], evtl. später vom Nutzer auswählabr ob berechneter Wert oder experimentell bestimmter Wert

    fct = 3.03

    def getFct(self):
        return self.fct

    #Ec = 22000 * (
    #        0.1 * fcm) ** 0.3  # mittlerer E-Modul aus Betondruckfestigkeit berechnet [N/mm²], evtl später vom Nutzer auswählbar ob berechneter Wert oder experimentell bestimmter Wert

    Ec = 22734

    def getEc(self):
        return self.Ec

    Es = 198800

    def getEs(self):
        return self.Es

    fy = 575  # Streckgrenze der Bewehrung [N/mm²]

    def getFy(self):
        return self.fy

    epsilony = fy / Es  # Dehnung bei Erreichen der Streckgrenze [-]

    def getEpsilony(self):
        return self.epsilony

    fu = 664  # Zugfestigkeit der Bewehrung [N/mm²]

    def getFu(self):
        return self.fu

    epsilonu = 0.1  # Bruchdehnung [-]

    def getEpsilonu(self):
        return self.epsilonu

    dag = 16  # Größtkorndurchmesser [mm]

    def getDag(self):
        return self.dag

    rhol = 0.01047  # Längsbewehrungsgrad [-] Später für Anwender: Auswahl ob Bestimmung der Fläche der Bewehrung über Längsbewehrungsgrad oder Angabe der Stäbe getroffen wird

    def getRhol(self):
        return self.rhol

    ds = 16  # Stabdurchmesser [mm]

    def getDs(self):
        return self.ds

    n = 3  # Stabanzahl [-]

    def getN(self):
        return self.n

    As = n * pi * (ds ** 2) / 4

    def getAs(self):
        return self.As

    alpha = 0.3   # Verhältniswert des oberen zum unteren Rissastes

    def getAlpha(self):
        return self.alpha

    # Eingabe von Schätzwerten für die Parameter deltaY, beta, sigmaz0 und phi
    # Alle Angaben sind in [mm], [°], [N/mm²] anzugeben
    # Eingabe erfolgt zurzeit in bruteforce_local.py daher hier auskommentiert!
    #deltaY = 1  # Rissfortschritt in vertikaler Richtung - MUSS kleiner werden [mm]

    #def getDeltaY(self):
    #    return self.deltaY

    #beta = radians(89)

    #def getBeta(self):
    #    return self.beta

    #sigmaz0 = 0.1  # Druckkraft aus Zahnbiegung - MUSS größer werden [N/mm²]

    #def getSigmaz0(self):
    #    return self.sigmaz0

    #phi = 0.0001  # Stabverdrehungswinkel - MUSS größer werden [°]

    #def getPhi(self):
    #    return self.phi

    # Eingabe der Startwerte von x0, beta, beta1, beta2, y1, y2 und y
    # [mm] [°]
    # x0 = d * 0.48  # Startwert der Betondruckzone zu Beginn der Iteration (Steuerungsgröße)

    #def getX0(self):
    #    return self.x0

    y1vor = 0.3  # Startwert des vertikalen Längenanteils des Abschnitt Bs des Risses, Unkonformität zu alpha = 0.3?

    def gety1vor(self):
        return self.y1vor

    y2vor = 0.7  # Startwert des vertikalen Längenanteils des Abschnitt As des Risses, Unkonformität zu alpha = 0.3?

    def gety2vor(self):
        return self.y2vor

    yvor = y1vor + y2vor  # Startwert der vertikalen Länge des Risses

    def getyvor(self):
        return self.yvor

    betavor = radians(90)

    def getBetavor(self):
        return self.betavor

    beta1vor = radians(90)  # Startwert der Rissneigung Abschnitt B, evtl. zu 90° annehmen?

    def getBeta1vor(self):
        return self.beta1vor

    beta2vor = radians(90)  # Startwert der Rissneigung Abschnitt A, evtl zu 90° annehmen?

    def getBeta2vor(self):
        return self.beta2vor

    # Startwerte h1vor und h2vor werden für Berechnungen nicht benötigt
