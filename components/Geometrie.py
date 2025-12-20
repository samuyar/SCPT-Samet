from math import *


class Geometrie:

    def __init__(self, d, yVor, y1Vor, y2Vor, deltaY, alpha, betavor, beta1vor, beta2vor):
        """
        Initiierung der Klasse Geometrie, erstellt/berechnet die Parameter des Rissfortschrittes
        :param d           statische Nutzhöhe
        :param yVor        Vertikale Risslänge aus vorherigem Schritt
        :param y1Vor       Vertikale Risslänge des Abschnitt B
        :param y2Vor       Vertikale Risslänge des Abschnitt A
        :param deltaY      Risszuwachs
        :param alpha       = 0.3, Vereinfachung der Anteile von y1 und y2 von y
        :param beta        Aktuell angesetzte Neigung des Risses an der Rissspitze (Zur Vereinfachten Abrufung des aktuellen beta-Wertes in den darauffolgenden Iterationen)
        :param betaVor     Neigung des Risss an der Rissspitze aus vorherigem Schritt
        :param beta1vor    Rissneigung Abschnitt B
        :param beta2vor    Rissneigung Abschnitt A
        Anmerkungen: Definitionen evtl. im späteren Verlauf in Getter/Setter auslagern?
        """
        self.scr = 0.7 * d  # Breite des Betonzahns nach Classen

        self.y = yVor + deltaY  # Neue vertikale Gesamtrisslänge

        self.y1 = alpha * self.y  # Vertikale Risslänge Abschnitt B
        self.y2 = (1 - alpha) * self.y  # Vertikale Risslänge Abschnitt A

        self.betavor = betavor  # Aktuell angesetzte Rissneigung an der Rissspitze

        # Entscheidungsfunktion für die Bestimmung der horizontalen Längenanteile
        self.h2 = None  # horizontale Risslänge des Abschnitt A
        if yVor > self.y2:
            self.h2 = 1 / tan(beta1vor) * y1Vor + 1 / tan(beta2vor) * y2Vor - (yVor - self.y2) / tan(beta1vor)
        else:
            self.h2 = 1 / tan(beta1vor) * y1Vor + 1 / tan(beta2vor) * y2Vor - (yVor - self.y2) / tan(betavor)

        self.h1 = None  # horizontale Risslänge Abschnit B
        if yVor > self.y2:
            self.h1 = 1 / tan(beta1vor) * y1Vor + 1 / tan(beta2vor) * y2Vor - self.h2 + deltaY / tan(betavor)
        else:
            self.h1 = self.y1 / tan(betavor)

        self.beta1 = atan(self.y1 / self.h1)  # Rissneigung Abschnitt B
        self.beta2 = atan(self.y2 / self.h2)  # Rissneigung Abschnitt A

        self.h3 = 1 / tan(self.beta1) * self.y1 + 1 / tan(self.beta2) * self.y2  # Länge des Delaminieriungsrisses

