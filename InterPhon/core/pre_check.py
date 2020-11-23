import numpy as np
from typing import List, Dict
from InterPhon import error


def to_int_numpy(data):
    try:
        if isinstance(data, int):
            data = [data, ]
        elif isinstance(data, str):
            data = [int(val) for val in data.strip().split()]
        else:
            data = [int(val) for val in data]
    except ValueError:
        raise ValueError("The items of '{0}' cannot be converted to int".format(data))
    return np.array(data)


class PreArgument(object):
    """
    Pre argument class to construct argument object during pre-process.
    The information about user arguments is stored in instance variables which are determined by
    'set_user_argument' method and validity is checked by 'check_user_argument' method.
    """
    def __init__(self, displacement: float = None,
                 enlargement: np.ndarray = None,
                 periodicity: np.ndarray = None):
        """
        Constructor of PreArgument class.

        :param displacement: (float) Displacement length (unit: Angst).
        :param enlargement: (np.ndarray[int]) Extension ratio along each a, b, c lattice.
        :param periodicity: (np.ndarray[bool]) Periodic (True) or not (False) along each a, b, c direction.
        """
        self.__displacement = displacement
        self.__enlargement = enlargement
        self.__periodicity = periodicity

    @property
    def displacement(self):
        return self.__displacement

    @displacement.setter
    def displacement(self, _displacement):
        _displacement = float(_displacement)

        if _displacement < 0.0:
            ValueError("Displacement length (unit: Angst) should be positive")
        elif _displacement > 0.1:
            print("Caution: the appropriate displacement (unit: Angst) is between 0.01 and 0.10")
            self.__displacement = _displacement
        else:
            self.__displacement = _displacement

    @property
    def enlargement(self):
        return self.__enlargement

    @enlargement.setter
    def enlargement(self, _enlargement):
        _enlargement = to_int_numpy(_enlargement)

        if len(_enlargement) != 3:
            raise error.Insufficient_DIM_Error
        else:
            self.__enlargement = _enlargement

    @property
    def periodicity(self):
        return self.__periodicity

    @periodicity.setter
    def periodicity(self, _periodicity):
        periodicity = []
        for value in _periodicity.strip().split():
            periodicity.append(True if value[0] not in ('0', 'F', 'f') else False)
        _periodicity = to_int_numpy(periodicity)

        if len(_periodicity) != 3:
            raise error.Insufficient_PBC_Error
        else:
            for ind, value in enumerate(_periodicity):
                if not value:
                    if self.enlargement[ind] != 1:
                        raise error.Mismatch_DIM_and_PBC_Error
            self.__periodicity = _periodicity

    def initialization(self):
        self.__displacement = None
        self.__enlargement = None
        self.__periodicity = None

    def set_user_argument(self, dict_args: Dict) -> None:
        """
        Method of PreArgument class.
        Set the variables of PreArgument instance from the information given by user.

        usage:
        " >>> instance_of_PreArgument.set_user_argument(list_args=arguments)"

        :param dict_args: (Dict) Arguments given by user.
        :return: (None)
        """
        self.initialization()
        for key, value in dict_args.items():
            if 'displacement' in key:
                self.displacement = value
            elif 'enlargement' in key:
                self.enlargement = value
            elif 'periodicity' in key:
                self.periodicity = value

    def check_user_argument(self) -> None:
        """
        Method of PreArgument class.
        Check the validity of instance variable.

        usage:
        " >>> instance_of_PreArgument.check_user_argument()"

        :return: (None)
        """
        if len(self.enlargement) != 3:
            raise error.Insufficient_DIM_Error

        if len(self.periodicity) != 3:
            raise error.Insufficient_PBC_Error

        for ind, value in enumerate(self.periodicity):
            if not value:
                if self.enlargement[ind] != 1:
                    raise error.Mismatch_DIM_and_PBC_Error

    def __deepcopy__(self, memodict: dict = {}) -> object:
        import copy

        cls = self.__class__
        result = cls.__new__(cls)
        memodict[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, copy.deepcopy(value, memodict))
        return result
