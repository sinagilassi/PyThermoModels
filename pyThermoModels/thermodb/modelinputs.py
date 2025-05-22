# ThermoModels
# import libs
import numpy as np
from typing import List, Dict, Union
from pyThermoDB import TableMatrixData
# local


class InputController:
    """
    InputController class for managing the inputs of the thermodynamic models including NRTL, UNIQUAC.
    """

    def __init__(self):
        """
        Initialize the InputController class.
        """
        pass

    def inputs_NRTL(self,
                    components: List[str],
                    activity_inputs: Dict[
                        str,
                        Union[List[float],
                              List[List[float]],
                              np.ndarray,
                              TableMatrixData]
                    ],
                    **kwargs):
        '''
        NRTL activity model for calculating activity coefficients.

        Parameters
        ----------
        components : list
            List of component names.
        activity_inputs : dict
            Dictionary containing the activity model inputs.
            - dg_ij : list, optional
                Interaction energy parameters for the components.
            - a_ij : list, optional
                Interaction energy parameters for the components.
            - b_ij : list, optional
                Interaction energy parameters for the components.
            - c_ij : list, optional
                Interaction energy parameters for the components.
            - d_ij : list, optional
                Interaction energy parameters for the components.
            - alpha_ij : list, optional
                Non-randomness parameter for the components.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.

        Returns
        -------
        inputs : dict
            Dictionary containing the inputs for the NRTL activity model.
            - dg_ij : list
                Interaction energy parameters for the components.
            - a_ij : list
                Interaction energy parameters for the components.
            - b_ij : list
                Interaction energy parameters for the components.
            - c_ij : list
                Interaction energy parameters for the components.
            - d_ij : list
                Interaction energy parameters for the components.
            - alpha_ij : list
                Non-randomness parameter for the components.
        '''
        try:
            # SECTION: analyze inputs
            # check
            if activity_inputs is None:
                raise ValueError(
                    "activity_inputs must be provided for NRTL activity model.")

            # NOTE: check if activity_inputs is a dictionary
            if activity_inputs is not None:
                # check if activity_inputs is a dictionary
                if not isinstance(activity_inputs, dict):
                    raise ValueError(
                        "activity_inputs must be a dictionary.")
                # check if activity_inputs is empty
                if len(activity_inputs) == 0:
                    raise ValueError(
                        "activity_inputs cannot be empty.")

            # NOTE: method 1
            # Δg_ij, interaction energy parameter
            dg_ij_src = activity_inputs.get(
                'dg_ij', None) or activity_inputs.get('dg', None)

            # NOTE: method 2
            # constants a, b, and c
            a_ij_src = activity_inputs.get(
                'a_ij', None) or activity_inputs.get('a', None)
            b_ij_src = activity_inputs.get(
                'b_ij', None) or activity_inputs.get('b', None)
            c_ij_src = activity_inputs.get(
                'c_ij', None) or activity_inputs.get('c', None)
            d_ij_src = activity_inputs.get(
                'd_ij', None) or activity_inputs.get('d', None)

            # NOTE: α_ij, non-randomness parameter
            alpha_ij_src = activity_inputs.get(
                'alpha_ij', None) or activity_inputs.get('alpha', None)
            if alpha_ij_src is None:
                raise ValueError(
                    "No valid source provided for non-randomness parameter (α_ij).")

            # NOTE: check method
            tau_ij_cal_method = 0

            # check
            if dg_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if a_ij_src is None or b_ij_src is None or c_ij_src is None or d_ij_src is None:
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dg_ij_src is not None:
                # use dg_ij
                if isinstance(dg_ij_src, TableMatrixData):
                    dg_ij = dg_ij_src.mat('dg', components)
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
            # α_ij, non-randomness parameter
            if isinstance(alpha_ij_src, TableMatrixData):
                alpha_ij = alpha_ij_src.mat('alpha', components)
            elif isinstance(alpha_ij_src, list):
                alpha_ij = np.array(alpha_ij_src)
            elif isinstance(alpha_ij_src, np.ndarray):
                alpha_ij = alpha_ij_src
            else:
                raise ValueError(
                    "Invalid source for non-randomness parameter (α_ij). Must be TableMatrixData, list of lists, or numpy array.")

            # SECTION: init inputs
            # set inputs
            inputs = {}

            # SECTION: set all binary interaction parameter matrix
            if tau_ij_cal_method == 1:
                # set
                inputs['dg_ij'] = dg_ij
            elif tau_ij_cal_method == 2:
                # set
                inputs['a_ij'] = a_ij
                inputs['b_ij'] = b_ij
                inputs['c_ij'] = c_ij
                inputs['d_ij'] = d_ij
            else:
                raise ValueError(
                    "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # SECTION: set α_ij
            # set
            inputs['alpha_ij'] = alpha_ij

            # res
            return inputs
        except Exception as e:
            raise Exception(f"Failed to generate NRTL inputs: {e}") from e

    def inputs_UNIQUAC(self,
                       components: List[str],
                       activity_inputs: Dict[
                           str,
                           Union[List[float],
                                 List[List[float]],
                                 np.ndarray,
                                 TableMatrixData]
                       ],
                       **kwargs):
        '''
        UNIQUAC activity model for calculating activity coefficients.

        Parameters
        ----------
        components : list
            List of component names.
        activity_inputs : dict
            Dictionary containing the activity model inputs.
            - dU_ij : list, optional
                Interaction energy parameters for the components.
            - a_ij : list, optional
                Interaction energy parameters for the components.
            - b_ij : list, optional
                Interaction energy parameters for the components.
            - c_ij : list, optional
                Interaction energy parameters for the components.
            - d_ij : list, optional
                Interaction energy parameters for the components.
            - r_i : list, optional
                Relative van der Waals volume of component i.
            - q_i : list, optional
                Relative van der Waals area of component i.
        kwargs : dict
            Additional parameters for the model.
            - interaction-energy-parameter : list, optional
                Interaction energy parameters for the components.

        Returns
        -------
        inputs : dict
            Dictionary containing the inputs for the UNIQUAC activity model.
            - dU_ij : list
                Interaction energy parameters for the components.
            - a_ij : list
                Interaction energy parameters for the components.
            - b_ij : list
                Interaction energy parameters for the components.
            - c_ij : list
                Interaction energy parameters for the components.
            - d_ij : list
                Interaction energy parameters for the components.
            - r_i : list
                Relative van der Waals volume of component i.
            - q_i : list
                Relative van der Waals area of component i.
        '''
        try:
            # SECTION: analyze inputs
            # check
            if activity_inputs is None:
                raise ValueError(
                    "activity_inputs must be provided for UNIQUAC activity model.")

            # NOTE: check if activity_inputs is a dictionary
            if activity_inputs is not None:
                # check if activity_inputs is a dictionary
                if not isinstance(activity_inputs, dict):
                    raise ValueError(
                        "activity_inputs must be a dictionary.")
                # check if activity_inputs is empty
                if len(activity_inputs) == 0:
                    raise ValueError(
                        "activity_inputs cannot be empty.")

            # NOTE: method 1
            # Δg_ij, interaction energy parameter
            dU_ij_src = activity_inputs.get(
                'dU_ij', None) or activity_inputs.get('dU', None)

            # NOTE: method 2
            # constants a, b, and c
            a_ij_src = activity_inputs.get(
                'a_ij', None) or activity_inputs.get('a', None)
            b_ij_src = activity_inputs.get(
                'b_ij', None) or activity_inputs.get('b', None)
            c_ij_src = activity_inputs.get(
                'c_ij', None) or activity_inputs.get('c', None)
            d_ij_src = activity_inputs.get(
                'd_ij', None) or activity_inputs.get('d', None)

            # NOTE: r_i, relative van der Waals volume of component i
            r_i_src = activity_inputs.get(
                'r_i', None) or activity_inputs.get('r', None)
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
            q_i_src = activity_inputs.get(
                'q_i', None) or activity_inputs.get('q', None)
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

            # NOTE: check method
            tau_ij_cal_method = 0
            if dU_ij_src is None:
                # check if a_ij, b_ij, c_ij are provided
                if (a_ij_src is None or
                    b_ij_src is None or
                    c_ij_src is None or
                        d_ij_src is None):
                    raise ValueError(
                        "No valid source provided for interaction energy parameter (Δg_ij) or constants a, b, c, and d.")
                # set method
                tau_ij_cal_method = 2

                # ! a_ij
                if isinstance(a_ij_src, TableMatrixData):
                    a_ij = a_ij_src.mat('a', components)
                elif isinstance(a_ij_src, list):
                    a_ij = np.array(a_ij_src)
                elif isinstance(a_ij_src, np.ndarray):
                    a_ij = a_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (a_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! b_ij
                if isinstance(b_ij_src, TableMatrixData):
                    b_ij = b_ij_src.mat('b', components)
                elif isinstance(b_ij_src, list):
                    b_ij = np.array(b_ij_src)
                elif isinstance(b_ij_src, np.ndarray):
                    b_ij = b_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (b_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! c_ij
                if isinstance(c_ij_src, TableMatrixData):
                    c_ij = c_ij_src.mat('c', components)
                elif isinstance(c_ij_src, list):
                    c_ij = np.array(c_ij_src)
                elif isinstance(c_ij_src, np.ndarray):
                    c_ij = c_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (c_ij). Must be TableMatrixData, list of lists, or numpy array.")

                # ! d_ij
                if isinstance(d_ij_src, TableMatrixData):
                    d_ij = d_ij_src.mat('d', components)
                elif isinstance(d_ij_src, list):
                    d_ij = np.array(d_ij_src)
                elif isinstance(d_ij_src, np.ndarray):
                    d_ij = d_ij_src
                else:
                    raise ValueError(
                        "Invalid source for interaction energy parameter (d_ij). Must be TableMatrixData, list of lists, or numpy array.")
            elif dU_ij_src is not None:
                # use dU_ij
                if isinstance(dU_ij_src, TableMatrixData):
                    dU_ij = dU_ij_src.mat('dU', components)
                elif isinstance(dU_ij_src, list):
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

            # SECTION: init inputs
            # set inputs
            inputs = {}

            # NOTE: set the binary interaction parameter matrix
            if tau_ij_cal_method == 1:
                # set
                inputs['dU_ij'] = dU_ij
            elif tau_ij_cal_method == 2:
                # set
                inputs['a_ij'] = a_ij
                inputs['b_ij'] = b_ij
                inputs['c_ij'] = c_ij
                inputs['d_ij'] = d_ij
            else:
                raise ValueError(
                    "Invalid method for calculating tau_ij. Must be 1 or 2.")

            # SECTION: set r_i and q_i
            # set
            inputs['r_i'] = r_i
            inputs['q_i'] = q_i

            # res
            return inputs
        except Exception as e:
            raise Exception(
                f"Failed to generate UNIQUAC inputs: {e}") from e
