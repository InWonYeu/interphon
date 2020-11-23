class Insufficient_PBC_Error(Exception):
    """
    Defined error class to guide for insufficient PBC arguments.
    """
    def __str__(self):
        """
        Error message for insufficient PBC arguments.

        :return: (str) Error message.
        """
        return "PBC is not sufficient, the length of PBC must be 3."


class Insufficient_DIM_Error(Exception):
    """
    Defined error class to guide for insufficient DIM arguments.
    """
    def __str__(self):
        """
        Error message for insufficient DIM arguments.

        :return: (str) Error message.
        """
        return "DIM is not sufficient, the length of DIM must be 3."


class Mismatch_DIM_and_PBC_Error(Exception):
    """
    Defined error class to guide for inconsistency between DIM and PBC arguments.
    """
    def __str__(self):
        """
        Error message for inconsistency between DIM and PBC arguments.

        :return: (str) Error message.
        """
        return "DIM is not matched with PBC. Super cell enlargement along the non-periodic direction is meaningless"


class Mismatch_DIM_post_Error(Exception):
    """
    Defined error class to guide for inconsistency between pre- and post- DIM arguments.
    """
    def __str__(self):
        """
        Error message for inconsistency between pre- and post- DIM arguments.

        :return: (str) Error message.
        """
        return "DIM of post-process is not matched with the DIM of pre-process."


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
    def __str__(self):
        """
        Error message for the attempt to access unspecified k-points.

        :return: (str) Error message.
        """
        return "The given k_point was not set."
