import numpy as np
from InterPhon import error
from .pre_check import PreArgument
from .unit_cell import UnitCell
from .super_cell import SuperCell


class PostArgument(PreArgument):
    """
    Post argument class to construct argument object during post-process.
    This child class is inherited from the PreArgument parent class.
    The information about user arguments is stored in instance variables which are determined by
    inherited 'set_user_argument' method.
    The validity of arguments is checked by the inherited 'check_user_argument' and the 'check_match_argument' method.
    """
    def __init__(self, displacement: float = None,
                 enlargement: np.ndarray = None,
                 periodicity: np.ndarray = None):
        """
        Constructor of PostArgument class.

        :param displacement: (float) Displacement length (unit: Angst).
        :param enlargement: (np.ndarray[int]) Extension ratio along each a, b, c lattice.
        :param periodicity: (np.ndarray[bool]) Periodic (True) or not (False) along each a, b, c direction.
        """
        super(PostArgument, self).__init__(displacement,
                                           enlargement,
                                           periodicity)

    def check_match_argument(self, unit_cell: UnitCell, super_cell: SuperCell) -> None:
        """
        Method of PostArgument class.
        Check the consistency of arguments between pre- and post-processes.

        usage:
        " >>> instance_of_PostArgument.check_match_argument(unit_cell=self.unit_cell, super_cell=self.super_cell)"

        :param unit_cell: (instance) of UnitCell class.
        :param super_cell: (instance) of SuperCell class.
        :return: (None)
        """
        for i in range(3):
            if not np.allclose(unit_cell.lattice_matrix[i, 0:3] * self.enlargement[i],
                               super_cell.lattice_matrix[i, 0:3]):
                raise error.Mismatch_DIM_post_Error

        _enlarge = 1
        for ind, value in enumerate(self.periodicity):
            if value:
                _enlarge = _enlarge * self.enlargement[ind]

        # for _unit_num, _super_num in zip(unit_cell.num_atom, super_cell.num_atom):
        #     if _unit_num * _enlarge != _super_num:
        #         raise error.Mismatch_DIM_post_Error
