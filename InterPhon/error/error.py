class Insufficient_PBC_Error(Exception):
    """
    An error class defined to announce insufficient PBC arguments.

    :param pbc: Periodicity
    :type pbc: :class:`core.PreArgument.periodicity`
    """
    def __init__(self, pbc):
        self.pbc = pbc

    def __str__(self):
        """
        Error message for insufficient PBC arguments.

        :return: Error message
        :rtype: str
        """
        return "The given PBC, {0}, is not sufficient. The length of PBC must be 3.".format(self.pbc)


class Insufficient_ENLARGE_Error(Exception):
    """
    An error class defined to announce insufficient ENLARGE arguments.

    :param enlarge: Enlargement
    :type enlarge: :class:`core.PreArgument.enlargement`
    """
    def __init__(self, enlarge):
        self.enlarge = enlarge

    def __str__(self):
        """
        Error message for insufficient ENLARGE arguments.

        :return: Error message
        :rtype: str
        """
        return "The given ENLARGE, {0}, is not sufficient. The length of ENLARGE must be 3.".format(self.enlarge)


class Mismatch_ENLARGE_and_PBC_Error(Exception):
    """
    An error class defined to announce inconsistency between ENLARGE and PBC arguments.

    :param enlarge: Enlargement
    :type enlarge: :class:`core.PreArgument.enlargement`
    :param pbc: Periodicity
    :type pbc: :class:`core.PreArgument.periodicity`
    """
    def __init__(self, enlarge, pbc):
        self.enlarge = enlarge
        self.pbc = pbc

    def __str__(self):
        """
        Error message for inconsistency between ENLARGE and PBC arguments.

        :return: Error message
        :rtype: str
        """
        return "ENLARGE, {0}, is not matched with PBC, {1}. " \
               "Supercell enlargement along the non-periodic direction is meaningless.".format(self.enlarge, self.pbc)


class Mismatch_ENLARGE_post_Error(Exception):
    """
    An error class defined to announce inconsistency between pre- and post- ENLARGE arguments.
    """
    def __str__(self):
        """
        Error message for inconsistency between pre- and post- ENLARGE arguments.

        :return: Error message
        :rtype: str
        """
        return "ENLARGE of post-process is not matched with the ENLARGE of pre-process."


class Mismatch_Kpath_and_PBC_Error(Exception):
    """
    An error class defined to announce inconsistency between K-path and PBC arguments.

    :param kpoint: K-point
    :type kpoint: np.ndarray[float]
    :param pbc: Periodicity
    :type pbc: :class:`core.PostArgument.periodicity`
    """
    def __init__(self, kpoint, pbc):
        self.kpoint = kpoint
        self.pbc = pbc

    def __str__(self):
        """
        Error message for inconsistency between K-path and PBC arguments.

        :return: Error message
        :rtype: str
        """
        return "K-point, {0}, is not matched with PBC, {1}. ".format(self.kpoint, self.pbc) + \
               "\nNon-zero reciprocal lattice along the non-periodic direction is meaningless."


class Invalid_Line_Kpath_Error(Exception):
    """
    An error class defined to announce invalid k-points setting for line path of band plot.
    """
    def __str__(self):
        """
        Error message for invalid k-points setting for line path of band plot.

        :return: Error message
        :rtype: str
        """
        return "The number of special k-points must be even."


class Not_Specified_Kpath_Error(Exception):
    """
    An error class defined to announce the attempt to access unspecified k-points.

    :param kpoint: K-point
    :type kpoint: np.ndarray[float]
    """
    def __init__(self, kpoint):
        self.kpoint = kpoint

    def __str__(self):
        """
        Error message for the attempt to access unspecified k-points.

        :return: Error message
        :rtype: str
        """
        return "The given k_point, {0}, was not set before.".format(self.kpoint)


class Cannot_Search_Point_Group(Exception):
    """
    An error class defined to announce that point group cannot be searched.

    :param look_up_table: Look-up-table used to identify point group (crystal class)
    :type look_up_table: np.ndarray[int]
    """
    def __init__(self, look_up_table):
        self.look_up_table = look_up_table

    def __str__(self):
        """
        Error message for the case that point group cannot be searched.

        :return: Error message
        :rtype: str
        """
        return "What is this point group? " \
               "Since the point group cannot be found, the symmetry function is turned off instead."


class Thermal_Imaginary_Frequency(Exception):
    """
    An error class defined to warn against an attempt to compute thermal properties with imaginary frequency.
    """
    def __str__(self):
        """
        Error message for an attempt to compute thermal properties with imaginary frequency.

        :return: Error message
        :rtype: str
        """
        return "Computation of thermal properties by imaginary frequency, which cannot be defined, is being attempted." \
               "\nThermal properties will be calculated by neglecting the corresponding imaginary frequency, " \
               "\nBut please be careful when using the results. "
