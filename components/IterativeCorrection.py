from components.ScptInput import CrackGeom
from stepOutput import StepOutput
from ScptVariableSet import ScptVariableSet
from cmath import atan, pi, sin, sqrt, tan
from math import cos
from scipy import optimize

class IterativeCorrection(StepOutput):
    """
    Contains output from equilibrium calculations
    """
    def __init__(self,
    phi: float,
    sigmaZ0: float,
    deltaY: float,
    beta: float
    ) -> None:
        super().__init__()
        self.phi = phi
        self.sigmaZ0 = sigmaZ0
        self.deltaY = deltaY
        self.beta = beta


def calculate(input: ScptVariableSet):
    """
    Calculations of initially estimated values for estimators
    - compatible reinforcement strain and rotation (eq. (61)-(66))
    - vertical stress at crack tip (eq. (67)-(70))
    - rotation of the member in the shear crack (eq. (71)-(73))
    - vertical crack expansion (eq. (74))
    """
    #Create object of class IterativeCorrection to store calculation results
    iterOutput = IterativeCorrection()


    #Effect of tension stiffening and rebar delamination: Determination of compatible reinforcement strain and rotation
    #Tensile force in reinforcement (eq. (61)) and tensile strain in long. reinf. (eq. (62) or (63)):
    Fs = input.equiOutput.Me / input.constitutiveOutput.z + input.constitutiveOutput.FaiPa * cos(input.geomOutput.beta2) + input.constitutiveOutput.FaiOr * (sin(input.geomOutput.beta2) - input.constitutiveOutput.zai / input.constitutiveOutput.z) + input.constitutiveOutput.Ffpz * (input.constitutiveOutput.zfpz / input.constitutiveOutput.z - sin(input.geomOutput.beta1)) + input.constitutiveOutput.Vuncr * input.constitutiveOutput.zuncr / input.constitutiveOutput.z + input.loading.Ne * (1 - input.geomdat.d / (2 * input.constitutiveOutput.z))        #TODO: Check whether this yields the same results as Fs from horizontal equilibrium

    #Reinforcement strain in control section is averaged value comprising different strain contributions
    #Contribution of pure reinforcement steel along delamination crack epsDelam:
    if (Fs / input.geomdat.As) < input.matpropReinf.fy:
        epsDelam = Fs / (input.matpropReinf.Es * input.geomdat.As)
    elif ((Fs / input.geomdat.As) >= input.matpropReinf.fy) ^ ((Fs / input.geomdat.As) < input.matpropReinf.fu):
        epsDelam = input.matpropReinf.fy / input.matpropReinf.Es + (Fs / input.geomdat.As - input.matpropReinf.fy) / ((input.matpropReinf.fu - input.matpropReinf.fy) / ((input.matpropReinf.fu - input.matpropReinf.fy) / input.matpropReinf.Es))
    else:
        ValueError('Strain in longitudinal reinforcement out of range.')

    #Contribution epsTs from regions with intact bond controlled by tension stiffening effect of surrounding concrete
    #Proposal by Bentz (eq. (64)):
    hceff = min(2, 5 * input.geomdat.d1, (input.geomdat.h - input.crack.x0) / 3)
    Aceff = hceff * input.geomdat.b
    #das ist falsch!!!!  M = Aceff / (pi * ds * n) wäre richtig statt:
    M = Aceff / input.geomdat.As

    def tsZero(eps: float) -> float:        #Here eps = epsTs
        eps * input.matpropReinf.Es + input.matpropConc.fct / (1 + sqrt(3.6 * M * eps)) - Fs / input.geomdat.As #eq. (64) (Fs/As must be equal to the right hand side of eq. (64))

    solution = optimize.newton(tsZero,0)    #Use Newton-Raphson iteration to solve for epsTs
    epsTs = solution[0]          #First entry in solution array is the estimated location where the function is zero

    #Calculation of averaged value of reinforcement strain epsS (eq. (65)):
    scrDelam = input.geomOutput.h3
    scr = input.geomOutput.scr

    if scrDelam <= scr:
        epsS = (epsTs * (scr - scrDelam) + epsDelam * scrDelam) / scr #eq. (65)
    else:
        epsS = epsDelam         #If delamination crack is present all along scr, then epsS = epsDelam (no bond)

    #Rotation of the member in the shear crack (eq. (66)):
    iterOutput.phi = epsS * scr / (input.geomdat.d - input.crack.x0) #eq. (66)


    #Cantilever action and bond of reinforcement: Determination of vertical stress at crack tip
    #Variation of tensile force between crack n and crack (n+1) (eq. (68)+(69):
    n = input.crack.lambdaCS * input.geomdat.d / scr
    deltaFsLim = (scr - scrDelam) * pi * (input.geomdat.ds * input.geomdat.n) * 1.26 * sqrt(input.matpropConc.fcm / 20) #bond strength according to ModelCode 2010 for ribbed bars #units are right (fcm in MPa and ds in mm)
    deltaFs = min((Fs * (n+1) / n - Fs), deltaFsLim) # -Fs is correct (mistake in the paper)

    #Clamping moment (eq. (67)):
    Mca = - input.constitutiveOutput.Ffpz * scr * cos(input.geomOutput.beta1) + input.constitutiveOutput.FaiOr * scr * cos(input.geomOutput.beta2) - input.constitutiveOutput.FaiPa * scr * sin(input.geomOutput.beta2) - input.constitutiveOutput.Vda * scr + deltaFs * (input.geomOutput.y2 + input.geomOutput.y1)

    #Vertical stress at the crack tip (eq. (70)):
    Wc = input.geomdat.b * scr ** 2 / 6                     #Assumption for moment of resistance of concrete tooth [mm³]
    iterOutput.sigmaZ0 = Mca / Wc


    #Crack propagation angle and (vertical) crack expansion (eq. (59), (60), ):
    #If initial estimation for beta was correct, calculated beta value should be the same. Otherwise: New iteration required.
    #Differentiation between shear span types:
    match input.structSys.shearSpanType:
        case 1:
            #Auxiliary equations for solution of quadratic equation for beta:
            u = 2 / 3 * input.geomdat.b * ((input.crack.x0 - input.crackPropOutput.x1) + 1 / 2 * ((input.crackPropOutput.x1 ** 2 + input.crack.x0 * input.crackPropOutput.x1) / (input.crack.x0 - input.crackPropOutput.x1))) * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0) * (input.constitutiveOutput.zuncr + input.crack.lambdaCS * input.geomdat.d)
            v = input.constitutiveOutput.Fcc * input.constitutiveOutput.z + input.constitutiveOutput.FaiOr * (input.constitutiveOutput.zai + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta2)) - input.constitutiveOutput.FaiPa * input.crack.lambdaCS * input.geomdat.d * sin(input.geomOutput.beta2) - input.constitutiveOutput.Ffpz * (input.constitutiveOutput.zfpz + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta1)) + input.loading.Ne * input.geomdat.d / 2- input.constitutiveOutput.Vda * input.crack.lambdaCS * input.geomdat.d - 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.crackPropOutput.sigma1 * input.constitutiveOutput.zct
            w = 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.constitutiveOutput.zct * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0)

            iterOutput.beta = [atan(u / (2 * v) + sqrt((u / (2 * v)) ** 2 - w / v)),
                               atan(u / (2 * v) - sqrt((u / (2 * v)) ** 2 - w / v))]
            iterOutput.beta.sort(reverse = True)
            iterOutput.beta = iterOutput.beta[0]    #greater value of beta is chosen as calculated beta #TODO: Check why greater value for beta is chosen

        case 2:
            #Auxiliary equations for solution of quadratic equation for beta: #TODO: Update with expression for sigmaBiaxFrac, x1 and zct depending on beta (numerical solution required, e.q. Newton-Raphson)
            u = 2 / 3 * input.geomdat.b * ((input.crack.x0 - input.crackPropOutput.x1) + 1 / 2 * ((input.crackPropOutput.x1 ** 2 + input.crack.x0 * input.crackPropOutput.x1) / (input.crack.x0 - input.crackPropOutput.x1))) * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0) * (input.constitutiveOutput.zuncr + input.crack.lambdaCS * input.geomdat.d * (1 + 1 / 2 * (input.crack.lambdaCS * input.geomdat.d) / (input.loading.a - input.crack.lambdaCS * input.geomdat.d)))
            v = input.constitutiveOutput.Fcc * input.constitutiveOutput.z + input.constitutiveOutput.FaiOr * (input.constitutiveOutput.zai + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta2) * (1 + 1 / 2 * (input.crack.lambdaCS * input.geomdat.d) / (input.loading.a - input.crack.lambdaCS * input.geomdat.d))) - input.constitutiveOutput.FaiPa * input.crack.lambdaCS * input.geomdat.d * sin(input.geomOutput.beta2) * (1 + 1 / 2 * (input.crack.lambdaCS * input.geomdat.d) / (input.loading.a - input.crack.lambdaCS * input.geomdat.d)) - input.constitutiveOutput.Ffpz * (input.constitutiveOutput.zfpz + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta1) * (1 + 1 / 2 * (input.crack.lambdaCS * input.geomdat.d) / (input.loading.a - input.crack.lambdaCS * input.geomdat.d))) + input.loading.Ne * input.geomdat.d / 2- input.constitutiveOutput.Vda * input.crack.lambdaCS * input.geomdat.d * (1 + 1 / 2 * (input.crack.lambdaCS * input.geomdat.d) / (input.loading.a - input.crack.lambdaCS * input.geomdat.d)) - 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.crackPropOutput.sigma1 * input.constitutiveOutput.zct
            w = 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.constitutiveOutput.zct * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0)

            iterOutput.beta = [atan(u / (2 * v) + sqrt((u / (2 * v)) ** 2 - w / v)),
                               atan(u / (2 * v) - sqrt((u / (2 * v)) ** 2 - w / v))]
            iterOutput.beta.sort(reverse = True)
            iterOutput.beta = iterOutput.beta[0]    #greater value of beta is chosen as calculated beta

        case 3:
            #Auxiliary equations for solution of quadratic equation for beta: #TODO: Update with expression for sigmaBiaxFrac, x1 and zct depending on beta (numerical solution required, e.q. Newton-Raphson)
            u = 2 / 3 * input.geomdat.b * ((input.crack.x0 - input.crackPropOutput.x1) + 1 / 2 * ((input.crackPropOutput.x1 ** 2 + input.crack.x0 * input.crackPropOutput.x1) / (input.crack.x0 - input.crackPropOutput.x1))) * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0) * (input.constitutiveOutput.zuncr + input.crack.lambdaCS * input.geomdat.d)
            v = input.constitutiveOutput.Fcc * input.constitutiveOutput.z + input.constitutiveOutput.FaiOr * (input.constitutiveOutput.zai + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta2)) - input.constitutiveOutput.FaiPa * input.crack.lambdaCS * input.geomdat.d * sin(input.geomOutput.beta2) - input.constitutiveOutput.Ffpz * (input.constitutiveOutput.zfpz + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta1)) + input.loading.Ne * input.geomdat.d / 2- input.constitutiveOutput.Vda * input.crack.lambdaCS * input.geomdat.d - 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.crackPropOutput.sigma1 * input.constitutiveOutput.zct - (input.structSys.a ** 2 + (input.crack.lambdaCS * input.geomdat.d) ** 2) / 2 * input.equiOutput.q
            w = 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.constitutiveOutput.zct * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0)

            iterOutput.beta = [atan(u / (2 * v) + sqrt((u / (2 * v)) ** 2 - w / v)),
                               atan(u / (2 * v) - sqrt((u / (2 * v)) ** 2 - w / v))]
            iterOutput.beta.sort(reverse = True)
            iterOutput.beta = iterOutput.beta[0]    #greater value of beta is chosen as calculated beta

        case 4:
            #Auxiliary equations for solution of quadratic equation for beta: #TODO: Update with expression for sigmaBiaxFrac, x1 and zct depending on beta (numerical solution required, e.q. Newton-Raphson)
            u = 2 / 3 * input.geomdat.b * ((input.crack.x0 - input.crackPropOutput.x1) + 1 / 2 * ((input.crackPropOutput.x1 ** 2 + input.crack.x0 * input.crackPropOutput.x1) / (input.crack.x0 - input.crackPropOutput.x1))) * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0) * (input.constitutiveOutput.zuncr + input.crack.lambdaCS * input.geomdat.d)
            v = input.constitutiveOutput.Fcc * input.constitutiveOutput.z + input.constitutiveOutput.FaiOr * (input.constitutiveOutput.zai + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta2)) - input.constitutiveOutput.FaiPa * input.crack.lambdaCS * input.geomdat.d * sin(input.geomOutput.beta2) - input.constitutiveOutput.Ffpz * (input.constitutiveOutput.zfpz + input.crack.lambdaCS * input.geomdat.d * cos(input.geomOutput.beta1)) + input.loading.Ne * input.geomdat.d / 2- input.constitutiveOutput.Vda * input.crack.lambdaCS * input.geomdat.d - 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.crackPropOutput.sigma1 * input.constitutiveOutput.zct - (input.structSys.a ** 2 + (input.crack.lambdaCS * input.geomdat.d) ** 2) / 2 * input.equiOutput.q - input.structSys.a * input.equiOutput.F
            w = 1 / 2 * input.geomdat.b * input.crackPropOutput.x1 * input.constitutiveOutput.zct * (input.crackPropOutput.sigma1 - iterOutput.sigmaZ0)

            iterOutput.beta = [atan(u / (2 * v) + sqrt((u / (2 * v)) ** 2 - w / v)),
                               atan(u / (2 * v) - sqrt((u / (2 * v)) ** 2 - w / v))]
            iterOutput.beta.sort(reverse = True)
            iterOutput.beta = iterOutput.beta[0]    #greater value of beta is chosen as calculated beta

        case _:
            raise ValueError("Chosen Shear Span Type not defined in this code.")



    #Vertical crack expansion (eq. (74)):
    iterOutput.deltaY = input.geomdat.d - input.crack.x0 - input.crack.yPrev - input.crackPropOutput.x1

    #Write iterative output to ScptVariableSet:
    input.iterOutput = iterOutput