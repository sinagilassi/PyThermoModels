# equation of states
# -------------------

# import packages/modules
import numpy as np
# internals
from ..configs import constants as CONST
from .eos import EOS


class EOSCore(EOS):
    # number of components
    componentsNo = 0
    # init

    def __init__(self, datasource, components, eos_name, mole_fraction, operating_conditions):
        self.datasource = datasource
        self.components = components
        self.eos_name = eos_name
        self.mole_fraction = mole_fraction
        params = operating_conditions
        # set PVT
        self.P = params.get("pressure", None)
        self.T = params.get("temperature", None)
        # check
        if self.P is None or self.T is None:
            raise Exception("P and T are required!")
        # component no
        self.componentsNo = len(self.components)
        # parent
        super().__init__(self.P, self.T, eos_name, mole_fraction)

    def _eosPR(self):
        '''
        Find compressibility factor (Z) at specified P and T then molar-volume (V).

        Parameters
        ----------
        None

        Returns
        -------
        res: dict
            Z: compressibility coefficient [-]
            eos-params: a,b,A,B,alpha,beta,gamma
        '''
        try:
            # component data
            componentsData = self.datasource

            # sorted data
            # ! Pc [Pa], Tc [K], w [-]
            # ! MPa to Pa
            componentsDataSorted = [
                [float(values['Pc']['value'])*1e6, float(
                    values['Tc']['value']), float(values['AcFa']['value'])]
                for item, values in componentsData.items()]

            # set a b matrix
            a = np.zeros(self.componentsNo)
            b = np.zeros(self.componentsNo)

            count = 0
            for item in componentsDataSorted:
                aLoop = self.aPR(item[0], item[1], item[2])
                # a.append(aLoop)
                a[count] = aLoop
                bLoop = self.bPR(item[0], item[1])
                # b.append(bLoop)
                b[count] = bLoop

                # T/Tc ratio
                T_Tc_ratio = self.T/item[1]

                # set
                count += 1

            # check pure, multi-component system
            if self.componentsNo > 1:
                # mixing rule to calculate a/b
                # kij
                kij = self.kijFill()
                aij = self.aijFill(a, kij)
                aSet = self.aMixing(aij, self.moleFraction)
                bSet = self.bMixing(b, self.moleFraction)
            else:
                # no change a/b (pure component)
                aSet = a[0]
                bSet = b[0]

            # set parameters A,B
            A = self.eos_A(aSet)
            B = self.eos_B(bSet)

            # build polynomial eos equation f(Z)
            alpha = self.eos_alpha(B)
            beta = self.eos_beta(A, B)
            gamma = self.eos_gamma(A, B)

            # eso-params
            esoParams = {
                "a": aSet,
                "b": bSet,
                "A": A,
                "B": B,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma
            }

            # find f(Z) root
            rootList = np.sort(self.findRootfZ(alpha, beta, gamma))

            # ! check how many real Z
            ZsNo = len(rootList)

            # z
            minZ = np.amin(rootList)
            maxZ = np.amax(rootList)

            # # molar-volume [m3/mol]
            # # -> all
            # molarVolumes = np.sort(self.molarVolume(rootList))
            # # -> liquid
            # molarVolumeLiquid = self.molarVolume(minZ)
            # # -> gas
            # molarVolumeGas = self.molarVolume(maxZ)

            # REVIEW
            # calculate specific volume [m^3/kg]
            # res
            res = {
                "Zs": rootList,
                "eos-params": esoParams,
                "pressure": self.P,
                "temperature": self.T
            }

            return res

        except Exception as e:
            raise Exception(e)

    def _eosMixPR(self):
        '''
        find compressibility factor (Z) at specified P and T
        then molar-volume is found.

        output:
            res:
                pressure: fixed pressure [Pa]
                temperature: fixed temperature [K]
                molar-volumes: for all Z [m^3/mol]
                gas: for the highest Z [m^3/mol]
                liquid: for the lowest Z [m^3/mol]
                Z: compressibility coefficient [-]
                eos-params: a,b,A,B,alpha,beta,gamma
        '''
        try:
            # component data
            componentsData = self.compData

            # sorted data
            # Pc [bar], Tc [K], w [-]
            # ! Pc [bar] => [Pa]
            componentsDataSorted = [
                [float(item['Pc'])*1e5, float(item['Tc']), float(item['w'])] for item in componentsData]

            # set a b matrix
            a = np.zeros(self.componentsNo)
            b = np.zeros(self.componentsNo)

            count = 0
            for item in componentsDataSorted:
                aLoop = self.aPR(item[0], item[1], item[2])
                # a.append(aLoop)
                a[count] = aLoop
                bLoop = self.bPR(item[0], item[1])
                # b.append(bLoop)
                b[count] = bLoop
                count += 1

            # check pure, multi-component system
            if self.componentsNo > 1:
                # mixing rule to calculate a/b
                # kij
                kij = self.kijFill()
                aij = self.aijFill(a, kij)
                aSet = self.aMixing(aij, self.moleFraction)
                bSet = self.bMixing(b, self.moleFraction)
            else:
                # no change a/b (pure component)
                aSet = a[0]
                bSet = b[0]

            # set parameters A,B
            A = self.eos_A(aSet)
            B = self.eos_B(bSet)

            # build polynomial eos equation f(Z)
            alpha = self.eos_alpha(B)
            beta = self.eos_beta(A, B)
            gamma = self.eos_gamma(A, B)

            # eso-params
            esoParams = {
                "a": aSet,
                "b": bSet,
                "A": A,
                "B": B,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma
            }

            # find f(Z) root
            rootList = self.findRootfZ(alpha, beta, gamma)

            # z
            minZ = np.amin(rootList)
            maxZ = np.amax(rootList)

            # molar-volume [m3/mol]
            # -> all
            molarVolumes = self.molarVolume(rootList)
            # -> liquid
            molarVolumeLiquid = self.molarVolume(minZ)
            # -> gas
            molarVolumeGas = self.molarVolume(maxZ)

            # REVIEW
            # calculate specific volume [m^3/kg]
            # res
            res = {
                "pressure": self.P,
                "temperature": self.T,
                "molar-volumes": molarVolumes*np.array(self.moleFraction),
                "gas": molarVolumeGas*np.array(self.moleFraction),
                "liquid": molarVolumeLiquid*np.array(self.moleFraction),
                "Z": rootList,
                "eos-params": esoParams
            }

            return res

        except Exception as e:
            raise Exception(e)

    def _eosVDW(self):
        '''
        van der Waals equation (VDW)

        find compressibility factor (Z) at specified P and T
        then molar-volume is found.

        output:
            res:
                Z: compressibility coefficient [-]
                eos-params: a,b,A,B,alpha,beta,gamma
        '''
        try:
            # component data
            componentsData = self.compData

            # sorted data
            # Pc [bar], Tc [K], w [-]
            # ! Pc [bar] => [Pa]
            componentsDataSorted = [
                [float(item['Pc'])*1e5, float(item['Tc']), float(item['w'])] for item in componentsData]

            # set a b matrix
            a = np.zeros(self.componentsNo)
            b = np.zeros(self.componentsNo)

            count = 0
            for item in componentsDataSorted:
                aLoop = self.aVDW(item[0], item[1])
                # a.append(aLoop)
                a[count] = aLoop
                bLoop = self.bVDW(item[0], item[1])
                # b.append(bLoop)
                b[count] = bLoop

                # T/Tc ratio
                T_Tc_ratio = self.T/item[1]

                # set
                count += 1

            # FIXME
            # check pure, multi-component system
            if self.componentsNo > 1:
                # mixing rule to calculate a/b
                # kij
                kij = self.kijFill()
                aij = self.aijFill(a, kij)
                aSet = self.aMixing(aij, self.moleFraction)
                bSet = self.bMixing(b, self.moleFraction)
            else:
                # no change a/b (pure component)
                aSet = a[0]
                bSet = b[0]

            # set parameters A,B
            A = self.eos_A(aSet)
            B = self.eos_B(bSet)

            # build polynomial eos equation f(Z)
            alpha = self.eos_alpha(B)
            beta = self.eos_beta(A, B)
            gamma = self.eos_gamma(A, B)

            # eso-params
            esoParams = {
                "a": aSet,
                "b": bSet,
                "A": A,
                "B": B,
                "alpha": alpha,
                "beta": beta,
                "gamma": gamma
            }

            # find f(Z) root
            rootList = np.sort(self.findRootfZ(alpha, beta, gamma))

            #! check how many real Z
            ZsNo = len(rootList)

            # z
            minZ = np.amin(rootList)
            maxZ = np.amax(rootList)

            # # molar-volume [m3/mol]
            # # -> all
            # molarVolumes = np.sort(self.molarVolume(rootList))
            # # -> liquid
            # molarVolumeLiquid = self.molarVolume(minZ)
            # # -> gas
            # molarVolumeGas = self.molarVolume(maxZ)

            # REVIEW
            # calculate specific volume [m^3/kg]
            # res
            res = {
                "Zs": rootList,
                "eos-params": esoParams,
                "pressure": self.P,
                "temperature": self.T
            }

            return res

        except Exception as e:
            raise Exception(e)
