from math import *
from scipy.integrate import quad


class Konstitutiv:

    def __init__(self, b, d, n, ds, Ec, fct, fcm, dag, phi, y1, y2, wfpz, r2, l2, tau0, sigma1, beta1, beta2,
                 sigmaX0, x0, x1, delta, deltaK, epsilonTop):


        # Ungerissene Druckzone
        self.Vuncr = 2 / 3 * b * tau0 * (
                    (x0 + x1) + 1 / 2 * ((x1) ** 2 + x0 * x1) / (x0 - x1))  # Querkraft in ungerissener Druckzone
        self.zuncr = y1 / tan(beta1) + y2 / tan(beta2)  # Hebelarm Vuncr

        # Zugkräfte
        self.Fct = 1 / 2 * b * sigmaX0 * x1  # Normalkraft (Zug)
        self.zct = y1 + y2 + 1 / 3 * x1  # Hebelarm Fct

      # Druckkraft der ungerissenen Betondruckzone
      #  ec1 = 2.2/1000
      #  k = 1.05 * Ec * (ec1 / fcm)

      #  def sigmac(shi, ec1, k, fcm):
      #      ec = shi * abs(epsilonTop)
      #      eta = ec / ec1
      #      return fcm * (k * eta - eta ** 2) / (1 + (k - 2) * eta)

      #  sigmac = quad(sigmac, 0, 1, args=(ec1, k, fcm))
      #  sigmac = sigmac[0]
      #  self.Fcc = abs(b * sigmac * x0)


        sigmaOK = max([epsilonTop * Ec, -fcm])
        self.Fcc = abs(
            1 / 2 * b * sigmaOK * x0)  # Druckkraft in ungerissener Druckzone, hier wurde in der Excel-Tabelle eine Fallunterscheidung in Abhängigkeit von x1 getroffen.

        # Bestimmung des Hebelarms von Fcc
        self.z = d - 1 / 3 * x0


        # Rissprozesszone
        #w1 = (0.028 * fcm ** 0.18 * dag ** 0.32) / sigma1
        w1 = (0.028 * fcm ** 0.18 * dag ** 0.32) / fct
        #self.Ffpz = b * y1 / sin(beta1) * sigma1 * w1 / wfpz * (1 - exp(- wfpz / w1))
        self.Ffpz = b * y1 / sin(beta1) * fct * w1 / wfpz * (1 - exp(- wfpz / w1))
        self.Vfpz = cos(beta1) * self.Ffpz
        zfpz1 = y1 / sin(beta1) * (1 - w1 / wfpz * (1 - (1 + wfpz / w1) * exp(- wfpz / w1)) / (1 - exp(- wfpz / w1)))
        self.zfpz = zfpz1 + y2 / sin(beta2) * cos(beta2 - beta1)

        self.wbot_0 = phi * r2
        self.wbot_1 = phi * (r2 + y2 / sin(beta2) * 1)

        # Rissreibung
        # Funktion für Rissschubspannungen (Aktualisiert 04/12/2021: Rough-Crack Model nach Gambarova & Karakoc (1983))
        def tauAi(shi, fcm, delta, phi, r2, y2, beta2, dag):
            w = phi * (r2 + y2 / sin(beta2) * shi)
            a_tau = 64 * (delta/dag)**1.3
            b_tau = 0.6/dag + (100 * w/dag) ** (2.4 + 75 * delta/dag)
            #r = delta / (phi * (r2 + y2 / sin(beta2) * shi))
            # a3 = 9.8 / fcm
            # a4 = 2.44 * (1 - 16 / fcm)
            # return max(0.25 * fcm * (1 - (2 * phi * (r2 + y2 / sin(beta2) * shi) / dag) ** 0.5) *
            #            r * (a3 + a4 * (abs(r)) ** 3) / (1 + a4 * r ** 4), 0)
            return (fcm ** 0.6) * a_tau / b_tau

        # Integration der Schubspannungen über die Risslänge
        tauAI = quad(tauAi, 0, 1, args=(fcm, delta, phi, r2, y2, beta2, dag))
        tauAI = tauAI[0]   # Erster Wert in Liste tauAI ist das korrekte Ergebnis für das Integral der Schubspannung
        self.FaiPa = b * l2 * tauAI  # Parallel zum Riss wirkende Kraft

        # Funktion der Rissspannungen (Aktualisiert 04/12/2021: Rough-Crack Model nach Gambarova & Karakoc (1983))
        def sigmaAi(shi, fcm, delta, phi, r2, y2, beta2, dag):
            w = phi * (r2 + y2 / sin(beta2) * shi)
            a_sigma = 2500 * (delta/dag)**2.3
            b_sigma = 0.6/dag + (100 * w/dag) ** (3 + 75 * delta/dag)
            # r = delta / (phi * (r2 + y2 / sin(beta2) * shi))
            # a3 = 9.8 / fcm
            # a4 = 2.44 * (1 - 16 / fcm)
            # return max(0.62 * (phi * (r2 + y2 / sin(beta2) * shi)) ** 0.5 * r / (1 + r ** 2) ** 0.25 *
            #            0.25 * fcm * (1 - (2 * phi * (r2 + y2 / sin(beta2) * shi) / dag) ** 0.5) *
            #            r * (a3 + a4 * (abs(r)) ** 3) / (1 + a4 * r ** 4), 0)
            return (fcm ** 0.6) * a_sigma / b_sigma

        # Integration der Spannungen über die Risslänge
        sigmaAI = quad(sigmaAi, 0, 1, args=(fcm, delta, phi, r2, y2, beta2, dag))
        sigmaAI = sigmaAI[0]  # Erster Wert in Liste sigmaAI ist das korrekte Ergebnis für das Integral der Spannungen
        self.FaiOr = sigmaAI * b * l2  # Orthogonal zum Riss wirkende Kraft

        # Bestimmen des Hebelarms zai von FaiOr durch Integration
        def zAi(shi, fcm, delta, phi, r2, y2, beta2, dag):
            # r = delta / (phi * (r2 + y2 / sin(beta2) * shi))
            # a3 = 9.8 / fcm
            # a4 = 2.44 * (1 - 16 / fcm)
            # return (0.62 * (phi * (r2 + y2 / sin(beta2) * shi)) ** 0.5 * r / (1 + r ** 2) ** 0.25 *
            #         0.25 * fcm * (1 - (2 * phi * (r2 + y2 / sin(beta2) * shi) / dag) ** 0.5) *
            #         r * (a3 + a4 * (abs(r)) ** 3) / (1 + a4 * r ** 4)) * (1 - shi)
            w = phi * (r2 + y2 / sin(beta2) * shi)
            a_sigma = 2500 * (delta / dag) ** 2.3
            b_sigma = 0.6 / dag + (100 * w / dag) ** (3 + 75 * delta / dag)
            return (fcm ** 0.6) * a_sigma / b_sigma * (1 - shi)

        zAI = quad(zAi, 0, 1, args=(fcm, delta, phi, r2, y2, beta2, dag))
        self.zai = None
        if self.FaiOr != 0:
            self.zai = (b * l2 ** 2) / self.FaiOr * zAI[0]
        else:
            self.zai = l2 / 2

        # Vertikalkomponente der Rissreibung
        self.Vai = self.FaiPa * sin(beta2) - self.FaiOr * cos(beta2)

        # Dübelwirkung
        bn = b - n * ds  # Breite des Betonquerschnitts auf Höhe der Bewehrung, Vereinfachung: Einlagige Bewehrung mit nur einem Stabdurchmesser, evtl. später schon in DefVar zu definieren und hier als konkreten Wert übergeben bekommen
        Vda0 = 1.64 * bn * ds * (fcm) ** (1 / 3)  # maximal aufnehmbare Querkraft durch Dübelwirkung
        # Entscheidungsfunktion zur Bestimmung der aufgenommenen Querkraft in Abhängigkeit von deltaK
        self.Vda = None  # resultierende aufgenommene Querkraft in Abhängigkeit von deltaK
        if deltaK < 0.05:
            self.Vda = Vda0 * (deltaK / 0.05) * (2 - deltaK / 0.05)
        elif deltaK > 0.05:
            self.Vda = max(Vda0 * (2.55 - deltaK) / 2.5, 0)

    def getVuncr(self):
        return self.Vuncr

    def getZuncr(self):
        return self.zuncr

    def getFct(self):
        return self.Fct

    def getZct(self):
        return self.zct

    def getFcc(self):
        return self.Fcc

    def getZ(self):
        return self.z

    def getFfpz(self):
        return self.Ffpz

    def getVfpz(self):
        return self.Vfpz

    def getZfpz(self):
        return self.zfpz

    def getFaiPa(self):
        return self.FaiPa

    def getFaiOr(self):
        return self.FaiOr

    def getVda(self):
        return self.Vda

    def getZai(self):
        return self.zai
