"""Projects module exceptions"""


class NotAProject(Exception):
    pass


class NotARun(Exception):
    pass


class NotDagError(Exception):
    """ Raised when an action would result in a network that is not a DAG """
    pass
