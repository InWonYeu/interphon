import numpy as np
from InterPhon import error
from .pre_check import PreArgument
from .unit_cell import UnitCell
from .super_cell import SuperCell


class PostArgument(PreArgument):
    """
    Post argument class to construct an argument object during post-process.
    This child class is inherited from the :class:`core.PreArgument` parent class.
    The information about user arguments is stored in the instance variables of this class.
    The instance variables are set by the inherited :class:`core.PreArgument.set_user_argument` method,
    and their validity is checked by the inherited :class:`core.PreArgument.check_user_argument` method and the :class:`core.PostArgument.check_match_argument` method.

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
        Constructor of PostArgument class.
        """
        super(PostArgument, self).__init__(displacement,
                                           enlargement,
                                           periodicity)

    def check_match_argument(self, unit_cell: UnitCell, super_cell: SuperCell) -> None:
        """
        Check the consistency of user arguments between pre- and post-processes.

        :param unit_cell: Instance of UnitCell class
        :type unit_cell: :class:`core.UnitCell`
        :param super_cell: Instance of SuperCell class
        :type super_cell: :class:`core.SuperCell`
        """
        for i in range(3):
            if not np.allclose(unit_cell.lattice_matrix[i, 0:3] * self.enlargement[i],
                               super_cell.lattice_matrix[i, 0:3]):
                raise error.Mismatch_ENLARGE_post_Error

        _enlarge = 1
        for ind, value in enumerate(self.periodicity):
            if value:
                _enlarge = _enlarge * self.enlargement[ind]

        # for _unit_num, _super_num in zip(unit_cell.num_atom, super_cell.num_atom):
        #     if _unit_num * _enlarge != _super_num:
        #         raise error.Mismatch_ENLARGE_post_Error
