from math import *
from scipy.integrate import quad


class Konstitutiv:

    def __init__(self, b, d, n, ds, Ec, fct, fcm, dag, phi, y1, y2, wfpz, r2, l2, tau0, sigma1, beta1, beta2,
                 sigmaX0, x0, x1, delta, deltaK, epsilonTop):

        # =========================
        # Parameter Grenzschichten
        # =========================
        p = 15.0  # Abstand der Grenzschichtmitten [mm]
        t = 2.0  # Grenzschichtdicke [mm]
        eta_inter = 0.95  # Reduktionsfaktor der Steifigkeit (E-Modul) in Interlayer #TEST wenn eta = 1.0, dann sollte es keine Unterschiede geben! (Ja, test passt!)

        def integrate_with_layers(stammfunktion, a0, a1):
            """
            Integriert s(x)*f(x) von a0 bis a1 über Grenzschichten.
            stammfunktion(x).
            s(x)=1 außerhalb, s(x)=eta innerhalb der Grenzschichten.
            """
            res = 0.0
            x_cur = a0

            # erste Grenzschichtmitte >= a0
            k = int(floor(a0 / p)) + 1

            while True:
                c = k * p
                a = c - t / 2.0
                bnd = c + t / 2.0
                if a >= a1:
                    break

                aa = max(a0, a)
                bb = min(a1, bnd)

                # normaler Bereich bis zur Grenzschicht
                if aa > x_cur:
                    res += (stammfunktion(aa) - stammfunktion(x_cur))

                # Grenzschichtbereich
                if bb > aa:
                    res += eta_inter * (stammfunktion(bb) - stammfunktion(aa))

                x_cur = max(x_cur, bb)
                k += 1

            # Rest normal
            if x_cur < a1:
                res += (stammfunktion(a1) - stammfunktion(x_cur))

            return res

        # =========================
        # Vuncr (Schubkraft, SCPT)
        # =========================
        L = x0 + x1
        denom = (x1 ** 2 - x0 ** 2)

        # Stammfunktion von tau_ref(x) = tau0/(x1^2-x0^2) * (x^2 - 2*x0*x)
        def F_tau(x):
            return (tau0 / denom) * (x ** 3 / 3.0 - x0 * x ** 2)

        self.Vuncr = b * integrate_with_layers(F_tau, 0.0, L)
        self.zuncr = y1 / tan(beta1) + y2 / tan(beta2)

        # =========================
        # Fct (Zugkraft, linear über x1 ab x0)
        # =========================
        def F_sig_t(x):
            return (sigmaX0 / (2.0 * x1)) * (x - x0) ** 2

        a = x0
        bnd = x0 + x1

        if bnd >= a:
            self.Fct = b * (F_sig_t(bnd) - F_sig_t(a))
        else:
            self.Fct = -b * (F_sig_t(a) - F_sig_t(bnd))

        self.zct = y1 + y2 + (1.0 / 3.0) * x1

        # =========================
        # Fcc (Druckkraft, linear über x0 ab Oberkante)
        # =========================
        sigmaOK = max([epsilonTop * Ec, -fcm])

        def F_sig_c(x):
            return sigmaOK * (x - (x ** 2) / (2.0 * x0))

        self.Fcc = abs(b * (F_sig_c(x0) - F_sig_c(0.0)))

        # Hebelarm
        self.z = d - (1.0 / 3.0) * x0


        # Rissprozesszone
        #w1 = (0.028 * fcm ** 0.18 * dag ** 0.32) / sigma1
        #w1 = (0.028 * fcm ** 0.18 * dag ** 0.32) / fct
        w1 = (0.04 * fcm ** 0.18 * dag ** 0.32) / fct # mit Gf für 3D-Beton

        #self.Ffpz = b * y1 / sin(beta1) * sigma1 * w1 / wfpz * (1 - exp(- wfpz / w1))
        #self.Ffpz = b * y1 / sin(beta1) * fct * w1 / wfpz * (1 - exp(- wfpz / w1))

        # Abminderungsfaktor nach Gl. (5-6): beta1 in [0, pi/2]
        beta_p = 1/1 # Abminderungswert für 3D-Druck, ergibt sich aus dem Verhältnis zwischen der Zugfestigkeit der Grenzfläche und der Matrix #TEST wenn = 1.0 dann sollte sich nichts ändern (Test bestanden)
        beta1_eff = min(max(beta1, 0.0), pi / 2.0)
        alpha_p_fpz = 1.0 - (4.0 * (1.0 - beta_p) / (pi ** 2)) * (beta1_eff - pi / 2.0) ** 2

        # Grundwert Ffpz
        Ffpz_ref = b * y1 / sin(beta1) * fct * w1 / wfpz * (1 - exp(- wfpz / w1))

        self.Ffpz = alpha_p_fpz * Ffpz_ref # Ffpz reduziert mit alpha_p_fpz

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
        cf = 1.0
        #self.FaiPa = b * l2 * tauAI  # Parallel zum Riss wirkende Kraft
        self.FaiPa = cf * b * l2 * tauAI  # Parallel zum Riss wirkende Kraft, mit Reduktionsfaktor

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
        #self.FaiOr = sigmaAI * b * l2  # Orthogonal zum Riss wirkende Kraft
        self.FaiOr = cf * sigmaAI * b * l2  # Orthogonal zum Riss wirkende Kraft, mit Reduktionsfaktor

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
        #Vda0 = 1.64 * bn * ds * (fcm) ** (1 / 3)  # maximal aufnehmbare Querkraft durch Dübelwirkung
        alpha_p_vda = 1 # Reduktion der Dübelwirkung für 3D-Druck
        Vda0 = 1.64 * bn * ds * (alpha_p_vda * fcm) ** (1 / 3)  # maximal aufnehmbare Querkraft durch Dübelwirkung

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
