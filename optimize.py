from components.DefVar import *
from components.Geometrie import *
from components.Gleichgewicht import *
from components.Kinematik import *
from components.Konstitutiv import *
from components.ParameterNeu import *
from components.Rissfortschritt import *
from math import *
from scipy.optimize import minimize, basinhopping


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

    def get_yvor_share(self):
        return self.yvor_share

    def get_y1vor_share(self):
        return self.y1vor_share

    def get_y2vor_share(self):
        return self.y2vor_share

    def get_betavor_share(self):
        return self.betavor_share

    def get_beta1vor_share(self):
        return self.beta1vor_share

    def get_beta2vor_share(self):
        return self.beta2vor_share

    def set_yvor_share(self, new_yvor):
        self.yvor_share = new_yvor

    def set_y1vor_share(self, new_y1vor):
        self.y1vor_share = new_y1vor

    def set_y2vor_share(self, new_y2vor):
        self.y2vor_share = new_y2vor

    def set_betavor_share(self, new_betavor):
        self.betavor_share = new_betavor

    def set_beta1vor_share(self, new_beta1vor):
        self.beta1vor_share = new_beta1vor

    def set_beta2vor_share(self, new_beta2vor):
        self.beta2vor_share = new_beta2vor


class SharingParamsEach:
    """
    Class for sharing parameters between calculation function an wrapper function for EACH run of optimization.
    """
    def __init__(self, beta, phi, sigma, deltaY, yvor, y1vor, y2vor, betavor, beta1vor, beta2vor):
        self.beta_share = beta
        self.phi_share = phi
        self.sigma_share = sigma
        self.deltaY_share = deltaY
        self.yvor_share = yvor
        self.y1vor_share = y1vor
        self.y2vor_share = y2vor
        self.betavor_share = betavor
        self.beta1vor_share = beta1vor
        self.beta2vor_share = beta2vor

    def get_beta_share(self):
        return self.beta_share

    def get_phi_share(self):
        return self.phi_share

    def get_sigma_share(self):
        return self.sigma_share

    def get_deltaY_share(self):
        return self.deltaY_share

    def get_yvor_share(self):
        return self.yvor_share

    def get_y1vor_share(self):
        return self.y1vor_share

    def get_y2vor_share(self):
        return self.y2vor_share

    def get_betavor_share(self):
        return self.betavor_share

    def get_beta1vor_share(self):
        return self.beta1vor_share

    def get_beta2vor_share(self):
        return self.beta2vor_share

    def set_beta_share(self, new_beta):
        self.beta_share = new_beta

    def set_phi_share(self, new_phi):
        self.phi_share = new_phi

    def set_sigma_share(self, new_sigma):
        self.sigma_share = new_sigma

    def set_deltaY_share(self, new_deltaY):
        self.deltaY_share = new_deltaY

    def set_yvor_share(self, new_yvor):
        self.yvor_share = new_yvor

    def set_y1vor_share(self, new_y1vor):
        self.y1vor_share = new_y1vor

    def set_y2vor_share(self, new_y2vor):
        self.y2vor_share = new_y2vor

    def set_betavor_share(self, new_betavor):
        self.betavor_share = new_betavor

    def set_beta1vor_share(self, new_beta1vor):
        self.beta1vor_share = new_beta1vor

    def set_beta2vor_share(self, new_beta2vor):
        self.beta2vor_share = new_beta2vor


def get_bounds():
    """
    Bounds for values.
    :return: List of bounds.
    """
    bounds_beta = (radians(80), radians(89))
    bounds_phi = (0.0001, 0.02)
    bounds_sigma = (-2.5, 2.5)
    bounds_delta = (1, 100)
    return bounds_beta, bounds_phi, bounds_sigma, bounds_delta


def calculate(opti_params, sharing_params, temp, x0, yvor, y1vor, y2vor, betavor, beta1vor, beta2vor):
    beta = opti_params[0]
    phi = opti_params[1]
    sigmaz0 = opti_params[2]
    deltaY = opti_params[3]

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
    ds = defVar.getDs()  # Stabdurchmesser [mm]
    n = defVar.getN()  # Anzahl der Stäbe [-]
    As = defVar.getAs()  # Fläche der Längsbewehrung [mm^2]
    alpha = defVar.getAlpha()  # Verhältniswert des oberen zum unteren Rissast
    Ne = defVar.getNe()  # Einwirkende Normalkraf
    # h1 = defVar.getH1vor()  # Vorherige Risslänge h1
    # h2 = defVar.getH2vor()  # Vorherige Risslänge h2

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
                           konstitutiv.getZct(), konstitutiv.getZuncr(), yvor)

        sumN = gleichgewicht.Fs - konstitutiv.Fcc + konstitutiv.Fct - konstitutiv.FaiOr * sin(
            geometrie.beta2) - konstitutiv.FaiPa * cos(geometrie.beta2) + konstitutiv.Ffpz * sin(
            geometrie.beta1) - Ne

        # Overwrite every run
        sharing_params.set_beta_share(neu.betaNeu)
        sharing_params.set_phi_share(neu.phiCalc)
        sharing_params.set_sigma_share(neu.sigmaZ0)
        sharing_params.set_deltaY_share(neu.deltaYneu)

        temp.set_betavor_share(neu.betaNeu)
        temp.set_beta1vor_share(geometrie.beta1)
        temp.set_beta2vor_share(geometrie.beta2)
        temp.set_yvor_share(geometrie.y)
        temp.set_y1vor_share(geometrie.y1)
        temp.set_y2vor_share(geometrie.y2)

        # print('{}, {}, {}, {}'.format(beta, phi, sigmaz0, deltaY))

    except:
        print('Error')

    return abs(sumN)


def get_init_values(sharing_params):
    return sharing_params.get_yvor_share(), sharing_params.get_y1vor_share(), sharing_params.get_y2vor_share(), sharing_params.get_beta2vor_share(), sharing_params.get_beta1vor_share(), sharing_params.get_beta2vor_share()


def run_optimizing():
    """
    Wrapper function
    :return: List of minimums.
    """
    def_var = DefVar()

    beta = radians(89)
    phi = 0.00015
    sigmaz0 = 0.1
    deltaY = 1
    initial_params = [beta, phi, sigmaz0, deltaY]

    sharing_params_each = SharingParamsEach(beta, phi, sigmaz0, deltaY, def_var.yvor, def_var.y1vor, def_var.y2vor, def_var.betavor, def_var.beta1vor, def_var.beta2vor)

    def rel_beta(initial_params):
        return (initial_params[0] / sharing_params_each.get_beta_share()) - 1

    def rel_phi(initial_params):
        return (initial_params[1] / sharing_params_each.get_phi_share()) - 1

    def rel_sigma(initial_params):
        return (initial_params[2] / sharing_params_each.get_sigma_share()) - 1

    def rel_deltaY(initial_params):
        return (initial_params[3] / sharing_params_each.get_deltaY_share()) - 1

    constraint_beta = {'type': 'eq', 'fun': rel_beta}
    constraint_phi = {'type': 'eq', 'fun': rel_phi}
    constraint_sigma = {'type': 'eq', 'fun': rel_sigma}
    constraint_deltaY = {'type': 'eq', 'fun': rel_deltaY}
    constraints = (constraint_beta, constraint_phi, constraint_sigma, constraint_deltaY)

    d = def_var.getD()
    x_list = [d * 0.5, d * 0.48, d * 0.46, d * 0.45, d * 0.43, d * 0.41, d * 0.39, d * 0.37, d * 0.35, d * 0.33,
                d * 0.31, d * 0.305, d * 0.30, d * 0.295, d * 0.29]

    #x_list = [d * 0.5]

    temp = Temp(def_var.yvor, def_var.y1vor, def_var.y2vor, def_var.betavor, def_var.beta1vor, def_var.beta2vor)

    for x0 in x_list:
        result = minimize(calculate, initial_params, args=(sharing_params_each, temp, x0, *get_init_values(sharing_params_each)), method='SLSQP', bounds=get_bounds(), tol=0.05, constraints=constraints)
        #result = basinhopping(calculate, initial_params, minimizer_kwargs={'method': 'SLSQP', 'args': (sharing_params_each, temp, x0, *get_init_values(sharing_params_each)), 'bounds': get_bounds(), 'tol': 0.05, 'constraints': constraints})

        print('Global Minimum: x = {}, f(x0) = {}'.format(result.x, result.fun))

        beta_neu = result['x'][0]
        phi_neu = result['x'][1]
        sigma_neu = result['x'][2]
        deltaY_neu = result['x'][3]
        initial_params = [beta_neu, phi_neu, sigma_neu, deltaY_neu]
        sharing_params_each.set_yvor_share(temp.get_yvor_share())
        sharing_params_each.set_y1vor_share(temp.get_y1vor_share())
        sharing_params_each.set_y2vor_share(temp.get_y2vor_share())
        sharing_params_each.set_betavor_share(temp.get_betavor_share())
        sharing_params_each.set_beta1vor_share(temp.get_beta1vor_share())
        sharing_params_each.set_beta2vor_share(temp.get_beta2vor_share())


run_optimizing()
