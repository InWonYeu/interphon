import numpy as np
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
    Pre argument class to construct an argument object during pre-process.
    The information about user arguments is stored in the instance variables of this class.
    The instance variables are set by the :class:`core.PreArgument.set_user_argument` method,
    and their validity is checked by the :class:`core.PreArgument.check_user_argument` method.

    :param displacement: Displacement length (unit: Angst), defaults to None
    :type displacement: float
    :param enlargement: Extension ratio along each a, b, c lattice direction, defaults to None
    :type enlargement: np.ndarray[int]
    :param periodicity: Periodic (True) or not (False) along each a, b, c direction, defaults to None
    :type periodicity: np.ndarray[bool]
    """
    def __init__(self, displacement: float = None,
                 enlargement: np.ndarray = None,
                 periodicity: np.ndarray = None):
        """
        Constructor of PreArgument class.
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
            raise error.Insufficient_ENLARGE_Error(_enlargement)
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
            raise error.Insufficient_PBC_Error(_periodicity)
        else:
            for ind, value in enumerate(_periodicity):
                if not value:
                    if self.enlargement[ind] != 1:
                        raise error.Mismatch_ENLARGE_and_PBC_Error(self.enlargement, _periodicity)
            self.__periodicity = _periodicity

    def initialization(self):
        """
        Initialize the instance variables.
        """
        self.__displacement = None
        self.__enlargement = None
        self.__periodicity = None

    def set_user_argument(self, dict_args: dict) -> None:
        """
        Set the variables of **PreArgument** instance from the information given by user.

        :param dict_args: Argument dictionary given by user
        :type dict_args: dict
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
        Check the validity of instance variable.
        """
        if len(self.enlargement) != 3:
            raise error.Insufficient_ENLARGE_Error(self.enlargement)

        if len(self.periodicity) != 3:
            raise error.Insufficient_PBC_Error(self.periodicity)

        for ind, value in enumerate(self.periodicity):
            if not value:
                if self.enlargement[ind] != 1:
                    raise error.Mismatch_ENLARGE_and_PBC_Error(self.enlargement, self.periodicity)

    def __deepcopy__(self, memodict: dict = {}) -> object:
        import copy

        cls = self.__class__
        result = cls.__new__(cls)
        memodict[id(self)] = result
        for key, value in self.__dict__.items():
            setattr(result, key, copy.deepcopy(value, memodict))
        return result
