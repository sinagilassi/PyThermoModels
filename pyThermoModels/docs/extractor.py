# ThermoModels
# import libs
import numpy as np
from typing import List, Dict, Optional
from pyThermoDB import TableMatrixData
from typing import Dict, List
import pycuc
# local


class Extractor:
    """
    Extractor class for managing thermodynamic models.
    """

    def __init__(
        self,
        components: List[str],
        datasource: Dict,
        equationsource: Dict
    ):
        """
        Initialize the Extractor class.

        Parameters
        ----------
        datasource: Dict
            Data source for the model
        equationsource: Dict
            Equation source for the model
        components: List[str]
            List of component names in the mixture

        """
        # SECTION: check parameters
        # Check datasource
        if not isinstance(datasource, dict):
            raise TypeError("datasource must be a dict")

        # Check equationsource
        if not isinstance(equationsource, dict):
            raise TypeError("equationsource must be a dict")

        # Check if components is a list
        if not isinstance(components, list):
            raise TypeError("components must be a list")

        # SECTION: Assign the parameters to instance variables
        self.datasource = datasource
        self.equationsource = equationsource
        self.components = components

    def NRTL_inputs(
        self,
        temperature: Optional[List[float | str]] = None,
        **kwargs
    ):
        '''
        Prepares inputs for the NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: check src
            # extract activity model inputs
            datasource = self.datasource

            # NOTE: check if datasource is a dictionary
            if datasource is not None:
                # check if datasource is a dictionary
                if not isinstance(datasource, dict):
                    raise ValueError(
                        "datasource must be a dictionary.")

                # check if datasource is empty
                if len(datasource) == 0:
                    raise ValueError(
                        "datasource cannot be empty.")

            # NOTE: check temperature
            if temperature is not None:
                # check if temperature is a list
                if not isinstance(temperature, list):
                    raise ValueError(
                        "temperature must be a list of floats or strings.")

                # set
                T_value = float(temperature[0])
                T_unit = str(temperature[1])

                # convert temperature to Kelvin
                T = pycuc.convert_from_to(
                    value=T_value,
                    from_unit=T_unit,
                    to_unit='K'
                )

            # NOTE: method 1
            # ! Δg_ij, interaction energy parameter
            dg_ij_src = datasource.get(
                'dg_ij', None) or datasource.get('dg', None)

            # NOTE: method 2
            # ! constants a, b, c, and d
            a_ij_src = datasource.get(
                'a_ij', None) or datasource.get('a', None)
            b_ij_src = datasource.get(
                'b_ij', None) or datasource.get('b', None)
            c_ij_src = datasource.get(
                'c_ij', None) or datasource.get('c', None)
            d_ij_src = datasource.get(
                'd_ij', None) or datasource.get('d', None)

            # NOTE: α_ij, non-randomness parameter
            alpha_ij_src = datasource.get(
                'alpha_ij', None) or datasource.get('alpha', None)
            if alpha_ij_src is None:
                # set default value
                alpha_ij_src = None

            # NOTE: tau_ij, binary interaction parameter
            tau_ij_src = datasource.get(
                'tau_ij', None) or datasource.get('tau', None)
            if tau_ij_src is None:
                # set default value
                tau_ij_src = None

            # SECTION: init NRTL model
            # activity model
            # activity_nrtl = ptm.activities(
            #     components=self.components, model_name='NRTL')

            # SECTION: extract data
            # NOTE: check method
            tau_ij_cal_method = 0
            if dg_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if a_ij_src is None or b_ij_src is None or c_ij_src is None or d_ij_src is None:
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', self.components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', self.components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', self.components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', self.components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dg_ij_src is not None:
                # ! use dg_ij
                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', self.components)
                elif isinstance(dg_ij_src, list):
                    dg_ij = np.array(dg_ij_src)
                elif isinstance(dg_ij_src, np.ndarray):
                    dg_ij = dg_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: extract data
            # NOTE: α_ij, non-randomness parameter
            # check
            if alpha_ij_src is not None or alpha_ij_src != 'None':
                if isinstance(alpha_ij_src, TableMatrixData):
                    alpha_ij = alpha_ij_src.mat('alpha', self.components)
                elif isinstance(alpha_ij_src, list):
                    alpha_ij = np.array(alpha_ij_src)
                elif isinstance(alpha_ij_src, np.ndarray):
                    alpha_ij = alpha_ij_src
                else:
                    alpha_ij = None
            else:
                # set default value
                alpha_ij = None

            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            # if tau_ij_cal_method == 1:
            #     tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M1(
            #         temperature=temperature, dg_ij=dg_ij)
            # elif tau_ij_cal_method == 2:
            #     tau_ij, tau_ij_comp = activity_nrtl.cal_tau_ij_M2(
            #         temperature=temperature,
            #         a_ij=a_ij, b_ij=b_ij, c_ij=c_ij, d_ij=d_ij)
            # else:
            #     raise ValueError(
            #         "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # NOTE: nrtl inputs
            inputs = {
                "alpha_ij": alpha_ij,
                "dg_ij": dg_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
                "temperature": temperature
            }

            # res
            return inputs
        except Exception as e:
            raise Exception(f"Failed to calculate NRTL activity: {e}") from e

    def UNIQUAC_inputs(
        self,
        temperature: Optional[List[float | str]] = None,
        **kwargs
    ):
        '''
        Prepares inputs for the UNIQUAC activity model for calculating activity coefficients.

        Parameters
        ----------
        temperature : List[float | str], optional
            Temperature in any units as: [300, 'K'], it is automatically converted to Kelvin.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.
        '''
        try:
            # SECTION: check src
            # extract activity model inputs
            datasource = kwargs.get('datasource', None)

            # NOTE: check if datasource is a dictionary
            if datasource is not None:
                # check if datasource is a dictionary
                if not isinstance(datasource, dict):
                    raise ValueError(
                        "datasource must be a dictionary.")
                # check if datasource is empty
                if len(datasource) == 0:
                    raise ValueError(
                        "datasource cannot be empty.")

            # NOTE: check temperature
            if temperature is not None:
                # check if temperature is a list
                if not isinstance(temperature, list):
                    raise ValueError(
                        "temperature must be a list of floats or strings.")

                # set
                T_value = float(temperature[0])
                T_unit = str(temperature[1])

                # convert temperature to Kelvin
                T = pycuc.convert_from_to(
                    value=T_value,
                    from_unit=T_unit,
                    to_unit='K'
                )

            # NOTE: method 1
            # Δg_ij, interaction energy parameter
            dU_ij_src = datasource.get(
                'dU_ij', None) or datasource.get('dU', None)

            # NOTE: method 2
            # constants a, b, and c
            a_ij_src = datasource.get(
                'a_ij', None) or datasource.get('a', None)
            b_ij_src = datasource.get(
                'b_ij', None) or datasource.get('b', None)
            c_ij_src = datasource.get(
                'c_ij', None) or datasource.get('c', None)
            d_ij_src = datasource.get(
                'd_ij', None) or datasource.get('d', None)

            # NOTE: r_i, relative van der Waals volume of component i
            r_i_src = datasource.get(
                'r_i', None) or datasource.get('r', None)
            if r_i_src is None:
                raise ValueError("No valid source provided for r_i.")

            # check if r_i is a list or numpy array
            if isinstance(r_i_src, list):
                r_i = np.array(r_i_src)
            elif isinstance(r_i_src, np.ndarray):
                r_i = r_i_src
            else:
                raise ValueError(
                    "Invalid source for r_i. Must be a list or numpy array.")

            # NOTE: q_i, relative van der Waals area of component i
            q_i_src = datasource.get(
                'q_i', None) or datasource.get('q', None)
            if q_i_src is None:
                raise ValueError("No valid source provided for q_i.")

            # check if q_i is a list or numpy array
            if isinstance(q_i_src, list):
                q_i = np.array(q_i_src)
            elif isinstance(q_i_src, np.ndarray):
                q_i = q_i_src
            else:
                raise ValueError(
                    "Invalid source for q_i. Must be a list or numpy array.")

            # SECTION: init NRTL model
            # activity model
            # activity_uniquac = ptm.activities(
            #     components=self.components, model_name='NRTL')

            # SECTION: extract data
            # NOTE: check method
            tau_ij_cal_method = 0
            if dU_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if a_ij_src is None or b_ij_src is None or c_ij_src is None or d_ij_src is None:
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', self.components)
                elif isinstance(a_ij_src, List[List[float]]):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', self.components)
                elif isinstance(b_ij_src, List[List[float]]):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', self.components)
                elif isinstance(c_ij_src, List[List[float]]):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', self.components)
                elif isinstance(d_ij_src, List[List[float]]):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dU_ij_src is not None:
                # ! use dU_ij
                if isinstance(dU_ij_src, TableMatrixData):
                    dU_ij = dU_ij_src.mat('dU', self.components)
                elif isinstance(dU_ij_src, List[List[float]]):
                    dU_ij = np.array(dU_ij_src)
                elif isinstance(dU_ij_src, np.ndarray):
                    dU_ij = dU_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (Δg_ij). Must be TableMatrixData, list of lists, or numpy array.")
                # set method
                tau_ij_cal_method = 1
            else:
                raise ValueError(
                    "No valid source provided for interaction energy parameter (Δg_ij) or constants A, B, C.")

            # SECTION: calculate tau_ij
            # NOTE: calculate the binary interaction parameter matrix (tau_ij)
            # if tau_ij_cal_method == 1:
            #     tau_ij, _ = activity_uniquac.cal_tau_ij_M1(
            #         temperature=temperature, dU_ij=dU_ij)
            # elif tau_ij_cal_method == 2:
            #     tau_ij, _ = activity_uniquac.cal_tau_ij_M2(
            #         temperature=temperature,
            #         a_ij=a_ij, b_ij=b_ij, c_ij=c_ij, d_ij=d_ij)
            # else:
            #     raise ValueError(
            #         "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # NOTE: nrtl inputs
            inputs = {
                "r_i": r_i,
                "q_i": q_i,
                "dU_ij": dU_ij,
                "a_ij": a_ij,
                "b_ij": b_ij,
                "c_ij": c_ij,
                "d_ij": d_ij,
                "temperature": temperature
            }

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to calculate UNIQUAC activity: {e}") from e
