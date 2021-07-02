import numpy as np
from typing import Union, Callable, Optional, Any, List

# Defined some data type
MatrixLike = List[Union[int, float]]
AtomType = List[str]
SelectIndex = List[int]
FilePath = Union[str, List[str]]
File = Union[None, List[None]]
KptPath = List[np.ndarray]
