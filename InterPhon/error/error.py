class Insufficient_PBC_Error(Exception):
    """
    Defined error class to guide for insufficient PBC arguments.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        """
        Error message for insufficient PBC arguments.

        :return: (str) Error message.
        """
        return "The given PBC, {0}, is not sufficient. The length of PBC must be 3.".format(self.value)


class Insufficient_ENLARGE_Error(Exception):
    """
    Defined error class to guide for insufficient ENLARGE arguments.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        """
        Error message for insufficient ENLARGE arguments.

        :return: (str) Error message.
        """
        return "The given ENLARGE, {0}, is not sufficient. The length of ENLARGE must be 3.".format(self.value)


class Mismatch_ENLARGE_and_PBC_Error(Exception):
    """
    Defined error class to guide for inconsistency between ENLARGE and PBC arguments.
    """
    def __init__(self, enlarge, pbc):
        self.enlarge = enlarge
        self.pbc = pbc

    def __str__(self):
        """
        Error message for inconsistency between ENLARGE and PBC arguments.

        :return: (str) Error message.
        """
        return "ENLARGE, {0}, is not matched with PBC, {1}. " \
               "Supercell enlargement along the non-periodic direction is meaningless.".format(self.enlarge, self.pbc)


class Mismatch_ENLARGE_post_Error(Exception):
    """
    Defined error class to guide for inconsistency between pre- and post- ENLARGE arguments.
    """
    def __str__(self):
        """
        Error message for inconsistency between pre- and post- ENLARGE arguments.

        :return: (str) Error message.
        """
        return "ENLARGE of post-process is not matched with the ENLARGE of pre-process."


class Invalid_Line_Kpath_Error(Exception):
    """
    Defined error class to guide for invalid k-points setting for line path of band plot.
    """
    def __str__(self):
        """
        Error message for invalid k-points setting for line path of band plot.

        :return: (str) Error message.
        """
        return "The number of special k-points must be even."


class Not_Specified_Kpath_Error(Exception):
    """
    Defined error class to guide for the attempt to access unspecified k-points.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        """
        Error message for the attempt to access unspecified k-points.

        :return: (str) Error message.
        """
        return "The given k_point, {0}, was not set before.".format(self.value)


class Cannot_Search_Point_Group(Exception):
    """
    Defined error class to notify that point group cannot be searched.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        """
        Error message for the case that point group cannot be searched.

        :return: (str) Error message.
        """
        return "What is this point group? " \
               "Since the point group cannot be found, the symmetry function is turned off instead."
