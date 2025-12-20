# Definition des Zugriffs auf eigene Pythondateien
from components.DefVar import *
from components.Geometrie import *
from components.Gleichgewicht import *
from components.Kinematik import *
from components.Konstitutiv import *
from components.ParameterNeu import *
from components.Rissfortschritt import *
from scipy.optimize import minimize
import math
# from scipy.optimize import root_scalar
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scipy.optimize import basinhopping


def calculate(variable):
    beta = variable[0]
    phi = variable[1]
    sigmaz0 = variable[2]
    deltaY = variable[3]


    defVar = DefVar()  # Erstellen des Objektes defVar zum Zugriff auf die Parameter (Hier dadurch verallgemeinert, spätere Erweiterung durch Benutzeroberfläche möglich)
    a = defVar.a  # halbe Spannweite des Systems [mm]
    F = defVar.F  # Einzellast auf statisches System [N]
    q = defVar.q  # Linienlast auf statisches System [N]
    d = defVar.getD()  # statische Nutzhöhe [mm]
    b = defVar.getB()  # Balkenbreite [mm]
    lambdaCS = defVar.lambdaCS  # Ort des Risses im Balken [mm]
    fcm = defVar.getFcm()  # mittlere Betondruckfestigkeit [N/mm^2]
    fct = defVar.getFct()  # Betonzugfestigkeit aus Betondruckfestigkeit berechnet [N/mm^2]
    Ec = defVar.getEc()  # E-Modul des Betons [N/mm^2]
    Es = defVar.getEs()  # E-Modul Stahl [N/mm^2]
    dag = defVar.getDag()  # Größtkorndurchmesser [mm]
    rhol = defVar.getRhol()  # Längsbewehrungsgrad [-]
    ds = defVar.getDs()  # Stabdurchmesser [mm]
    n = defVar.getN()  # Anzahl der Stäbe [-]
    As = defVar.getAs()  # Fläche der Längsbewehrung [mm^2]
    alpha = defVar.getAlpha()  # Verhältniswert des oberen zum unteren Rissast
    Ne = defVar.getNe()  # Einwirkende Normalkraft

    # Vorgegebene Startwerte zu Beginn der Iteration, Variablen werden zum Teil während der Iteration neu belegt
    x0 = defVar.getX0()  # Höhe der Betondruckzone zu Beginn der Iteration
    #phi = defVar.getPhi()  # Stabverdrehungswinkel - MUSS größer werden [Bogenmaß]
    #deltaY = defVar.getDeltaY()  # Rissfortschritt in vertikaler Richtung - MUSS kleiner werden [mm]
    #beta = defVar.getBeta()
    betavor = defVar.getBetavor()  # Rissfortschrittswinkel des vorherigen Schrittes beta(i-1)
    #sigmaz0 = defVar.getSigmaz0()  # Druckkraft aus Zahnbiegung - MUSS größer werden [N/mm²]
    y1vor = defVar.gety1vor()  # Vorheriger Wert des vertikalen Längenanteils des Abschnitt 1 des Risses
    y2vor = defVar.gety2vor()  # Vorheriger Wert des vertikalen Längenanteils des Abschnitt 2 des Risses
    yvor = defVar.getyvor()  # Vorheriger Wert der vertikalen Länge des Risses
    #h1 = defVar.getH1vor()  # Vorherige Risslänge h1
    #h2 = defVar.getH2vor()  # Vorherige Risslänge h2
    beta1vor = defVar.getBeta1vor()  # Vorheriger Wert der Rissneigung Abschnitt 1
    beta2vor = defVar.getBeta2vor()  # Vorheriger Wert der Rissneigung Abschnitt 2

    try:
        geometrie = Geometrie(d, yvor, y1vor, y2vor, deltaY, alpha, betavor, beta1vor, beta2vor)

        # Bestimmen des Rissfortschritts des Iterationsschritts
        rissfortschritt = Rissfortschritt(d, geometrie.scr, x0, geometrie.y1, geometrie.y2, beta, phi, fct, fcm,
                                          Ec,
                                          sigmaz0)

        # Bestimmen der Risskinematik
        kinematik = Kinematik(phi, geometrie.scr, geometrie.y1, geometrie.y2, x0, rissfortschritt.x1,
                              geometrie.beta1, geometrie.beta2)
        #        kinematik = Kinematik(phi, geometrie.scr, geometrie.y1, geometrie.y2, x0, rissfortschritt.getX1(), geometrie.beta1, geometrie.beta2)

        # Berechnung der Traganteile
        konstitutiv = Konstitutiv(b, d, n, ds, Ec, fct, fcm, dag, phi, geometrie.y1, geometrie.y2,
                                  kinematik.wfpz,
                                  kinematik.r2, kinematik.l2, rissfortschritt.getTau0(), geometrie.beta1,
                                  geometrie.beta2, rissfortschritt.sigmaX0, x0, rissfortschritt.getX1(),
                                  kinematik.getDelta(), kinematik.getDeltak(), kinematik.getEpsilonTop())

        # Bestimmen der resultierenden Stahlzugkraft, Querkraft, Normalkraft und des resultierenden Momentes
        gleichgewicht = Gleichgewicht(As, Es, a, d, q, F, konstitutiv.getFcc(), konstitutiv.getFct(),
                                      konstitutiv.getFaiPa(), konstitutiv.getFaiOr(), konstitutiv.getFfpz(),
                                      konstitutiv.getVuncr(), konstitutiv.getVda(), geometrie.beta1,
                                      geometrie.beta2, konstitutiv.getZ(), konstitutiv.getZct(),
                                      konstitutiv.getZai(), konstitutiv.getZfpz(), konstitutiv.getZuncr(),
                                      lambdaCS,
                                      rissfortschritt.getX1(), Ne)

        # Prüfen ob die Parameter-Annahmen für die Verdrehung und den Rissneigungswinkel korrekt waren, bzw erneute Bestimmung dieser
        neu = ParameterNeu(lambdaCS, geometrie.scr, d, b, phi, x0, rissfortschritt.getX1(), geometrie.y1,
                           geometrie.y2, deltaY, beta, geometrie.beta1, geometrie.beta2, gleichgewicht.epsilonS,
                           rissfortschritt.sigma1, Ne, gleichgewicht.Fs, konstitutiv.getFcc(),
                           konstitutiv.getFfpz(), konstitutiv.getFaiOr(), konstitutiv.getFaiPa(),
                           konstitutiv.getVda(), konstitutiv.getZ(), konstitutiv.getZai(),
                           konstitutiv.getZfpz(),
                           konstitutiv.getZct(), konstitutiv.getZuncr(), yvor, fcm)

        sumN = gleichgewicht.Fs - konstitutiv.Fcc + konstitutiv.Fct - konstitutiv.FaiOr * sin(
            geometrie.beta2) - konstitutiv.FaiPa * cos(geometrie.beta2) + konstitutiv.Ffpz * sin(
            geometrie.beta1) - Ne

        def rel_beta():
            return (variable[0] / neu.betaNeu) - 1

        def rel_phi():
            return (variable[1] / neu.phiCalc) - 1

        def rel_sigma():
            return (variable[2] / neu.sigmaZ0) - 1

        def rel_deltaY():
            return (variable[3] / neu.deltaYneu) - 1

        const_rel_beta = sqrt(rel_beta() ** 2)
        const_rel_phi = sqrt(rel_phi() ** 2)
        const_rel_sigma = sqrt(rel_sigma() ** 2)
        const_rel_deltaY = sqrt(rel_deltaY() ** 2)
        const_sumN = sqrt((sumN / konstitutiv.getFcc()) ** 2)

        ret_val = const_rel_beta + const_rel_phi + const_rel_sigma + const_rel_deltaY + const_sumN
        return ret_val

    except Exception as error:
        print('Math error')
        return 1000.


variable = [radians(89), 0.0001, 0, 1]

bounds_beta = (radians(15), radians(89))
bounds_phi = (0.0001, 0.02)
bounds_sigma = (-2.5, 2.5)
bounds_delta = (1, 100)

bnds = (bounds_beta, bounds_phi, bounds_sigma, bounds_delta)

sol = basinhopping(calculate, variable, minimizer_kwargs={"method": "SLSQP", "bounds": bnds}, stepsize=0.001)
print(sol)

