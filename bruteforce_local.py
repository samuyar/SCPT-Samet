from components.DefVar import *
from components.Geometrie import *
from components.Gleichgewicht import *
from components.Kinematik import *
from components.Konstitutiv import *
from components.ParameterNeu import *
from components.Rissfortschritt import *
from math import *
from scipy.optimize import brute

import pandas as pd
import matplotlib.pyplot as plt


class Temp:
    """
    Class for sharing parameters between calculation function an wrapper function for ONE run of optimization.
    """

    def __init__(self, yvor, y1vor, y2vor, betavor, beta1vor, beta2vor):
        self.yvor_share = yvor
        self.y1vor_share = y1vor
        self.y2vor_share = y2vor
        self.betavor_share = betavor
        self.beta1vor_share = beta1vor
        self.beta2vor_share = beta2vor

        # Zwischenspeichern von Werten für Plots
        self.x_beta = []
        self.y_relbeta = []
        self.x_phi = []
        self.y_relphi = []
        self.x_sigma = []
        self.y_relsigma = []
        self.x_deltaY = []
        self.y_reldeltaY = []
        #self.z_error = []

    def get_yvor_share(self):
        return self.yvor_share

    def get_y1vor_share(self):
        return self.y1vor_share

    def get_y2vor_share(self):
        return self.y2vor_share

    def get_h1vor_share(self):
        return self.h1vor_share

    def get_h2vor_share(self):
        return self.h2vor_share

    def get_h3vor_share(self):
        return self.h3vor_share

    def get_betavor_share(self):
        return self.betavor_share

    def get_beta1vor_share(self):
        return self.beta1vor_share

    def get_beta2vor_share(self):
        return self.beta2vor_share

    def get_x1_share(self):
        return self.x1_share

    def get_sumN_share(self):
        return self.sumN_share

    def get_Vuncr_share(self):
        return self.Vuncr_share

    def get_Vda_share(self):
        return self.Vda_share

    def get_Vai_share(self):
        return self.Vai_share

    def get_Vfpz_share(self):
        return self.Vfpz_share

    def get_Fcc_share(self):
        return self.Fcc_share

    def get_Fct_share(self):
        return self.Fct_share

    def get_FaiPa_share(self):
        return self.FaiPa_share

    def get_FaiOr_share(self):
        return self.FaiOr_share

    def get_Ffpz_share(self):
        return self.Ffpz_share

    def get_Fs_share(self):
        return self.Fs_share

    def get_err_share(self):
        return self.err_share

    def get_relbeta_share(self):
        return self.relbeta_share

    def get_relphi_share(self):
        return self.relphi_share

    def get_relsigma_share(self):
        return self.relsigma_share

    def get_reldeltay_share(self):
        return self.reldeltay_share

    def get_relN_share(self):
        return self.relN_share

    def get_Vges_share(self):
        return self.Vges_share

    def get_beta_calc_share(self):
        return self.beta_calc_share

    def get_phi_calc_share(self):
        return self.phi_calc_share

    def get_sigma_calc_share(self):
        return self.sigma_calc_share

    def get_deltaY_calc_share(self):
        return self.deltaY_calc_share

    def get_deltak_share(self):
        return self.deltak_share

    def get_deltabot_share(self):
        return self.deltabot_share

    def get_wbot1_share(self):
        return self.wbot1_share

    def get_wbot0_share(self):
        return self.wbot0_share

    def get_wfpz_share(self):
        return self.wfpz_share

    def get_epsilons_share(self):
        return self.epsilons_share

    def get_epsiloncr_share(self):
        return self.epsiloncr_share

    def get_epsilontop_share(self):
        return self.epsilontop_share

    def get_tau0_share(self):
        return self.tau0_share

    def get_sigmax0_share(self):
        return self.sigmax0_share

    def get_sigma1_share(self):
        return self.sigma1_share




    def set_yvor_share(self, new_yvor):
        self.yvor_share = new_yvor

    def set_y1vor_share(self, new_y1vor):
        self.y1vor_share = new_y1vor

    def set_y2vor_share(self, new_y2vor):
        self.y2vor_share = new_y2vor

    def set_h1vor_share(self, new_h1vor):
        self.h1vor_share = new_h1vor

    def set_h2vor_share(self, new_h2vor):
        self.h2vor_share = new_h2vor

    def set_h3vor_share(self, new_h3vor):
        self.h3vor_share = new_h3vor

    def set_betavor_share(self, new_betavor):
        self.betavor_share = new_betavor

    def set_beta1vor_share(self, new_beta1vor):
        self.beta1vor_share = new_beta1vor

    def set_beta2vor_share(self, new_beta2vor):
        self.beta2vor_share = new_beta2vor

    def set_x1_share(self, new_x1):
        self.x1_share = new_x1

    def set_sumN_share(self, new_sumN):
        self.sumN_share = new_sumN

    def set_Vuncr_share(self, new_Vuncr):
        self.Vuncr_share = new_Vuncr

    def set_Vda_share(self, new_Vda):
        self.Vda_share = new_Vda

    def set_Vai_share(self, new_Vai):
        self.Vai_share = new_Vai

    def set_Vfpz_share(self, new_Vfpz):
        self.Vfpz_share = new_Vfpz

    def set_Fcc_share(self, new_Fcc):
        self.Fcc_share = new_Fcc

    def set_Fct_share(self, new_Fct):
        self.Fct_share = new_Fct

    def set_FaiPa_share(self, new_FaiPa):
        self.FaiPa_share = new_FaiPa

    def set_FaiOr_share(self, new_FaiOr):
        self.FaiOr_share = new_FaiOr

    def set_Ffpz_share(self, new_Ffpz):
        self.Ffpz_share = new_Ffpz

    def set_Fs_share(self, new_Fs):
        self.Fs_share = new_Fs

    def set_err_share(self, new_err):
        self.err_share = new_err

    def set_relbeta_share(self, new_relbeta):
        self.relbeta_share = new_relbeta

    def set_relphi_share(self, new_relphi):
        self.relphi_share = new_relphi

    def set_relsigma_share(self, new_relsigma):
        self.relsigma_share = new_relsigma

    def set_reldeltay_share(self, new_reldeltay):
        self.reldeltay_share = new_reldeltay

    def set_relN_share(self, new_relN):
        self.relN_share = new_relN

    def set_Vges_share(self, new_Vges):
        self.Vges_share = new_Vges

    def set_beta_calc_share(self, new_beta_calc):
        self.beta_calc_share = new_beta_calc

    def set_phi_calc_share(self, new_phi_calc):
        self.phi_calc_share = new_phi_calc

    def set_sigma_calc_share(self, new_sigma_calc):
        self.sigma_calc_share = new_sigma_calc

    def set_deltaY_calc_share(self, new_deltaY_calc):
        self.deltaY_calc_share = new_deltaY_calc

    def set_deltak_share(self, new_deltak):
        self.deltak_share = new_deltak

    def set_deltabot_share(self, new_deltabot):
        self.deltabot_share = new_deltabot

    def set_wbot1_share(self, new_wbot1):
        self.wbot1_share = new_wbot1

    def set_wbot0_share(self, new_wbot0):
        self.wbot0_share = new_wbot0

    def set_wfpz_share(self, new_wfpz):
        self.wfpz_share = new_wfpz

    def set_epsilons_share(self, new_epsilons):
        self.epsilons_share = new_epsilons

    def set_epsiloncr_share(self, new_epsiloncr):
        self.epsiloncr_share = new_epsiloncr

    def set_epsilontop_share(self, new_epsilontop):
        self.epsilontop_share = new_epsilontop

    def set_tau0_share(self, new_tau0):
        self.tau0_share = new_tau0

    def set_sigmax0_share(self, new_sigmax0):
        self.sigmax0_share = new_sigmax0

    def set_sigma1_share(self, new_sigma1):
        self.sigma1_share = new_sigma1



def calculate(x, temp, x0, yvor, y1vor, y2vor, betavor, beta1vor, beta2vor):
    beta = x[0]
    phi = x[1]
    sigmaz0 = x[2]
    deltaY = x[3]

    defVar = DefVar()  # Erstellen des Objektes defVar zum Zugriff auf die Parameter (Hier dadurch verallgemeinert, spätere Erweiterung durch Benutzeroberfläche möglich)
    a = defVar.a  # halbe Spannweite des Systems [mm]
    F = defVar.F  # Einzellast auf statisches System [N]
    q = defVar.q  # Linienlast auf statisches System [N]
    d = defVar.getD()  # statische Nutzhöhe [mm]
    b = defVar.getB()  # Balkenbreite [mm]
    h = defVar.getH()  # Balkenhöhe [mm]
    lambdaCS = defVar.lambdaCS  # Ort des Risses im Balken [mm]
    fcm = defVar.getFcm()  # mittlere Betondruckfestigkeit [N/mm^2]
    fct = defVar.getFct()  # Betonzugfestigkeit aus Betondruckfestigkeit berechnet [N/mm^2]
    Ec = defVar.getEc()  # E-Modul des Betons [N/mm^2]
    Es = defVar.getEs()  # E-Modul Stahl [N/mm^2]
    fy = defVar.getFy()  # Streckgrenze [N/mm^2]
    epsilony = defVar.getEpsilony()  # Dehnung bei Streckgrenze [-]
    fu = defVar.getFu()  # Zugfestigkeit Stahl [-]
    epsilonu = defVar.getEpsilonu()  # Bruchdehnung Stahl [-]
    dag = defVar.getDag()  # Größtkorndurchmesser [mm]
    ds = defVar.getDs()  # Stabdurchmesser [mm]
    n = defVar.getN()  # Anzahl der Stäbe [-]
    As = defVar.getAs()  # Fläche der Längsbewehrung [mm^2]
    alpha = defVar.getAlpha()  # Verhältniswert des oberen zum unteren Rissast
    Ne = defVar.getNe()  # Einwirkende Normalkraf

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
                                  kinematik.r2, kinematik.l2, rissfortschritt.getTau0(), rissfortschritt.getSigma1(),
                                  geometrie.beta1, geometrie.beta2, rissfortschritt.sigmaX0, x0, rissfortschritt.getX1(),
                                  kinematik.getDelta(), kinematik.getDeltak(), kinematik.getEpsilonTop())

        # Bestimmen der resultierenden Stahlzugkraft, Querkraft, Normalkraft und des resultierenden Momentes
        gleichgewicht = Gleichgewicht(As, Es, a, d, q, F, konstitutiv.getFcc(), konstitutiv.getFct(),
                                      konstitutiv.getFaiPa(), konstitutiv.getFaiOr(), konstitutiv.getFfpz(),
                                      konstitutiv.getVuncr(), konstitutiv.getVda(), geometrie.beta1,
                                      geometrie.beta2, konstitutiv.getZ(), konstitutiv.getZct(),
                                      konstitutiv.getZai(), konstitutiv.getZfpz(), konstitutiv.getZuncr(),
                                      lambdaCS,
                                      rissfortschritt.getX1(), Ne, fy, fu, epsilony, epsilonu, geometrie.h3,
                                      geometrie.scr, b, h, x0, n, ds, fct, fcm)

        # Prüfen ob die Parameter-Annahmen für die Verdrehung und den Rissneigungswinkel korrekt waren, bzw erneute Bestimmung dieser
        neu = ParameterNeu(lambdaCS, geometrie.scr, d, b, h, phi, x0, rissfortschritt.getX1(), geometrie.y1,
                           geometrie.y2, deltaY, beta, geometrie.beta1, geometrie.beta2, gleichgewicht.epsilonS, gleichgewicht.epsilonSdelam,
                           rissfortschritt.sigma1, Ne, gleichgewicht.Fs, konstitutiv.getFcc(),
                           konstitutiv.getFfpz(), konstitutiv.getFaiOr(), konstitutiv.getFaiPa(),
                           konstitutiv.getVda(), konstitutiv.getZ(), konstitutiv.getZai(),
                           konstitutiv.getZfpz(),
                           konstitutiv.getZct(), konstitutiv.getZuncr(), yvor,
                           gleichgewicht.scrdelam, ds, n, fcm)

        sumN = gleichgewicht.Fs - konstitutiv.Fcc + konstitutiv.Fct - konstitutiv.FaiOr * sin(
            geometrie.beta2) - konstitutiv.FaiPa * cos(geometrie.beta2) + konstitutiv.Ffpz * sin(
            geometrie.beta1) - Ne

        relN = (sumN / konstitutiv.Fcc)

        try:
            rel_beta = neu.betaNeu / beta - 1
            rel_phi = neu.phiCalc / phi - 1
            rel_sigma = neu.sigmaZ0 / sigmaz0 - 1
            rel_deltaY = neu.deltaYneu / deltaY - 1

            err = (sqrt((rel_phi ** 2) + (rel_deltaY ** 2) + (rel_beta ** 2) + (rel_sigma ** 2)))

            # Übergabe der berechneten Werte an den temporären Speicher
            temp.set_betavor_share(neu.betaNeu)
            temp.set_beta1vor_share(geometrie.beta1)
            temp.set_beta2vor_share(geometrie.beta2)
            temp.set_yvor_share(geometrie.y)
            temp.set_y1vor_share(geometrie.y1)
            temp.set_y2vor_share(geometrie.y2)
            temp.set_h1vor_share(geometrie.h1)
            temp.set_h2vor_share(geometrie.h2)
            temp.set_h3vor_share(geometrie.h3)
            temp.set_x1_share(rissfortschritt.x1)
            temp.set_sumN_share(sumN)
            temp.set_Vuncr_share(konstitutiv.Vuncr)
            temp.set_Vda_share(konstitutiv.Vda)
            temp.set_Vai_share(konstitutiv.Vai)
            temp.set_Vfpz_share(konstitutiv.Vfpz)
            temp.set_Fcc_share(konstitutiv.Fcc)
            temp.set_Fct_share(konstitutiv.Fct)
            temp.set_FaiPa_share(konstitutiv.FaiPa)
            temp.set_FaiOr_share(konstitutiv.FaiOr)
            temp.set_Ffpz_share(konstitutiv.Ffpz)
            temp.set_Fs_share(gleichgewicht.Fs)
            temp.set_err_share(err)
            temp.set_relbeta_share(rel_beta)
            temp.set_relphi_share(rel_phi)
            temp.set_relsigma_share(rel_sigma)
            temp.set_reldeltay_share(rel_deltaY)
            temp.set_relN_share(relN)
            temp.set_Vges_share(gleichgewicht.Ve)
            temp.set_beta_calc_share(neu.betaNeu)
            temp.set_phi_calc_share(neu.phiCalc)
            temp.set_sigma_calc_share(neu.sigmaZ0)
            temp.set_deltaY_calc_share(neu.deltaYneu)
            temp.set_deltak_share(kinematik.deltak)
            temp.set_deltabot_share(kinematik.delta)
            temp.set_wbot1_share(konstitutiv.wbot_1)
            temp.set_wbot0_share(konstitutiv.wbot_0)
            temp.set_wfpz_share(kinematik.wfpz)
            temp.set_epsilons_share(gleichgewicht.epsilonS)
            temp.set_epsiloncr_share(kinematik.epsilonCr)
            temp.set_epsilontop_share(kinematik.epsilonTop)
            temp.set_tau0_share(rissfortschritt.tau0)
            temp.set_sigmax0_share(rissfortschritt.sigmaX0)
            temp.set_sigma1_share(rissfortschritt.sigma1)

            #print(err)
            #print(neu.sigmaZ0, neu.betaNeu, neu.deltaYneu, neu.phiCalc)
            #print(sigmaz0, degrees(beta), deltaY, phi)
            #print(sumN)
            #print(relN)
            #print(' ')

            # Speichern bzw. Übergabe von Werten für Plots
            temp.x_beta.append(beta)
            temp.y_relbeta.append(rel_beta)
            temp.x_phi.append(phi)
            temp.y_relphi.append(rel_phi)
            temp.x_sigma.append(sigmaz0)
            temp.y_relsigma.append(rel_sigma)
            temp.x_deltaY.append(deltaY)
            temp.y_reldeltaY.append(rel_deltaY)
            #temp.z_error.append(err)

        except AttributeError as error:
            err = 10
            print(error)

    except UnboundLocalError as error:
        print('Value error: {}'.format(error))
        print('Error')

    except OverflowError as error:
        print('Value error: {}'.format(error))
        print('Error')


    return err


def opti():

    index = 0
    result_df = pd.DataFrame(columns=['x0', 'Error', 'beta', 'phi', 'sigmaz0', 'deltaY',
                                      'y', 'y1', 'y2', 'h1', 'h2', 'h3', 'beta1', 'beta2', 'x1', 'sumN',
                                      'Vges', 'Vuncr', 'Vfpz', 'Vai', 'Vda',
                                      'Fcc', 'Fct', 'Ffpz', 'FaiPa', 'FaiOr', 'Fs',
                                      'rel_beta', 'rel_phi', 'rel_sigma', 'rel_deltaY', 'relN',
                                      'beta_calc', 'phi_calc', 'sigmaz0_calc', 'deltaY_calc',
                                      'wfpz','wbot0', 'wreinf', 'deltabot', 'deltak',
                                      'epsilons', 'epsiloncr', 'epsilontop', 'tau0', 'sigmax0', 'sigma1'])

    # Init Params
    def_var = DefVar()

    beta = radians(89)
    phi = 0.00005
    sigmaz0 = 0.02
    deltaY = 50

    initial_params = [beta, phi, sigmaz0, deltaY]

    temp = Temp(def_var.yvor, def_var.y1vor, def_var.y2vor, def_var.betavor, def_var.beta1vor, def_var.beta2vor)

    x_list = [def_var.d * 0.61,
              def_var.d * 0.605,
              def_var.d * 0.60,
              def_var.d * 0.595,
              def_var.d * 0.59,
              def_var.d * 0.585,
              def_var.d * 0.58,
              def_var.d * 0.575,
	          def_var.d * 0.57,
              def_var.d * 0.565,
              def_var.d * 0.56,
              def_var.d * 0.555,
              def_var.d * 0.55,
              def_var.d * 0.545,
	          def_var.d * 0.54,
              def_var.d * 0.535,
              def_var.d * 0.53,
              def_var.d * 0.525,
              def_var.d * 0.52,
              def_var.d * 0.515,
  	          def_var.d * 0.51,
	          def_var.d * 0.505,
	          def_var.d * 0.500,
              def_var.d * 0.495,
              def_var.d * 0.49]




    for x0_list in x_list:

        # beta
        ranges = (
        slice(max(initial_params[0] - radians(15), radians(0.1)), initial_params[0], radians(1)),
        # phi
        slice(max(initial_params[1] - 0.0001, 0.00005), initial_params[1] + 0.0006, 0.00001),
        # sigmaZ0
        slice(max(initial_params[2] - 0.05, 0.00), initial_params[2] + 0.50, 0.01),
        # deltaY
        slice(max(initial_params[3] - 50, 0.1), initial_params[3] + 20, 1))

        result = brute(calculate, ranges, args=(temp, x0_list,
                                                temp.get_yvor_share(),
                                                temp.get_y1vor_share(),
                                                temp.get_y2vor_share(),
                                                temp.get_betavor_share(),
                                                temp.get_beta1vor_share(),
                                                temp.get_beta2vor_share()),
                       full_output=False, workers=48)

        beta_neu = min(initial_params[0], result[0])
        phi_neu = max(initial_params[1], result[1])
        sigma_neu = result[2]
        deltaY_neu = max(result[3], 0)

        # Für Plots
        #fixed_sigma = max(initial_params[2] - 0.05, 0.020)
        #fixed_deltaY = max(initial_params[3] - 50, 1)

        initial_params = [beta_neu, phi_neu, sigma_neu, deltaY_neu]

        result_df.loc[index] = [x0_list, temp.get_err_share(),
                                beta_neu, phi_neu, sigma_neu, deltaY_neu,
                                temp.get_yvor_share(), temp.get_y1vor_share(), temp.get_y2vor_share(),
                                temp.get_h1vor_share(), temp.get_h2vor_share(), temp.get_h3vor_share(),
                                temp.get_beta1vor_share(), temp.get_beta2vor_share(),
                                temp.get_x1_share(), temp.get_sumN_share(), temp.get_Vges_share(),
                                temp.get_Vuncr_share(), temp.get_Vfpz_share(), temp.get_Vai_share(), temp.get_Vda_share(),
                                temp.get_Fcc_share(), temp.get_Fct_share(), temp.get_Ffpz_share(), temp.get_FaiPa_share(),
                                temp.get_FaiOr_share(), temp.get_Fs_share(),
                                temp.get_relbeta_share(), temp.get_relphi_share(), temp.get_relsigma_share(),
                                temp.get_reldeltay_share(), temp.get_relN_share(),
                                temp.get_beta_calc_share(), temp.get_phi_calc_share(), temp.get_sigma_calc_share(),
                                temp.get_deltaY_calc_share(), temp.get_wfpz_share(), temp.get_wbot0_share(),
                                temp.get_wbot1_share(), temp.get_deltabot_share(), temp.get_deltak_share(),
                                temp.get_epsilons_share(), temp.get_epsiloncr_share(), temp.get_epsilontop_share(),
                                temp.get_tau0_share(), temp.get_sigmax0_share(), temp.get_sigma1_share()]

        index += 1
        result_df.to_csv('results.csv', index=False)

        # Erstellen der Plots
        # Info: Plots werden nur richtig erzeugt, wenn mit einem Kern gerechnet wird (workers=1)
        # x_beta = temp.x_beta
        # y_relbeta = temp.y_relbeta
        # result_x_beta = x_beta[-1]
        # result_y_relbeta = y_relbeta[-1]
        #
        # x_phi = temp.x_phi
        # y_relphi = temp.y_relphi
        # result_x_phi = x_phi[-1]
        # result_y_relphi = y_relphi[-1]
        #
        # x_sigma = temp.x_sigma
        # y_relsigma = temp.y_relsigma
        # result_x_sigma = x_sigma[-1]
        # result_y_relsigma = y_relsigma[-1]
        #
        # x_deltaY = temp.x_deltaY
        # y_reldeltaY = temp.y_reldeltaY
        # result_x_deltaY = x_deltaY[-1]
        # result_y_reldeltaY = y_reldeltaY[-1]

        #z_error = temp.z_error
        #result_z_error = temp.z_error[-1]

        # fig_1 = plt.figure(1, figsize=(20, 10))
        # chart_1 = fig_1.add_subplot(221)
        # chart_2 = fig_1.add_subplot(222)
        # chart_3 = fig_1.add_subplot(223)
        # chart_4 = fig_1.add_subplot(224)
        #
        # chart_1.scatter(x_beta, y_relbeta, s=5)
        # chart_1.scatter(result_x_beta, result_y_relbeta, s=5)
        # chart_1.set_xlabel('beta schätz')
        # chart_1.set_ylabel('rel_beta')
        #
        # chart_2.scatter(x_phi, y_relphi, s=5)
        # chart_2.scatter(result_x_phi, result_y_relphi, s=5)
        # chart_2.set_xlabel('phi schätz')
        # chart_2.set_ylabel('rel_phi')
        #
        # chart_3.scatter(x_sigma, y_relsigma, s=5)
        # chart_3.scatter(result_x_sigma, result_y_relsigma, s=5)
        # chart_3.set_xlabel('sigma schätz')
        # chart_3.set_ylabel('rel_sigma')
        #
        # chart_4.scatter(x_deltaY, y_reldeltaY, s=5)
        # chart_4.scatter(result_x_deltaY, result_y_reldeltaY, s=5)
        # chart_4.set_xlabel('deltaY schätz')
        # chart_4.set_ylabel('rel_deltaY')
        #
        # plt.savefig("Grenzen" + str(x0_list) + ".png")

        #indices = []
        #extract_x_beta = []
        #extract_x_phi = []
        #extract_z_err = []

        #for index1, index2 in enumerate(x_sigma):
        #    if index2 == fixed_sigma and x_deltaY[index1] == fixed_deltaY:
        #        indices.append(index1)
        #        extract_x_beta.append(x_beta[index1])
        #        extract_x_phi.append(x_phi[index1])
        #        extract_z_err.append(z_error[index1])

        #fig_2 = plt.figure(2)
        #ax2 = fig_2.add_subplot(111, projection='3d')
        #ax2.scatter3D(extract_x_beta, extract_x_phi, extract_z_err, c=extract_z_err, cmap='viridis')
        #ax2.scatter3D(result_x_beta, result_x_phi, result_z_error, c='red')
        #ax2.set_xlabel('beta schätz')
        #ax2.set_ylabel('phi schätz')
        #ax2.set_zlabel('err')

        #fig_3 = plt.figure(3)
        #ax3 = fig_3.add_subplot(111, projection='3d')
        #ax3.scatter3D(x_beta, x_phi, z_error, c=z_error, cmap='viridis')
        #ax3.scatter3D(result_x_beta, result_x_phi, result_z_error, c='red')
        #ax3.set_xlabel('beta schätz')
        #ax3.set_ylabel('phi schätz')
        #ax3.set_zlabel('err')

        #plt.show()
        # fig_1.clear()
        #
        # del temp.x_beta[:]
        # del temp.y_relbeta[:]
        # del temp.x_phi[:]
        # del temp.y_relphi[:]
        # del temp.x_sigma[:]
        # del temp.y_relsigma[:]
        # del temp.x_deltaY[:]
        # del temp.y_reldeltaY[:]
        #del temp.z_error[:]


if __name__ == "__main__":
    opti()
