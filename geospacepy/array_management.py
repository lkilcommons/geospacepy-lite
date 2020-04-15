import numpy as np

def is_number_or_len_one_array(number):
    """Deterimine if a given value is either, one, a python numerical primative
    (int or float), or two, a numpy array of length one.

    PARAMETERS
    ----------

        number, float,int or np.ndarray
            The object to check for singluarness
    
    RETURNS
    -------

        is_single_value, bool
            True if the object is a float or int, or a numpy array of length
            one. False otherwise
    """
    try:
        l = len(number)
        if l != 1 and isinstance(number,np.ndarray):
            return False
        else:
            return True
    except TypeError:
        if isinstance(number,int) or isinstance(number,float):
            return True
        else:
            return False

def check_follows_3_component_vector_convention(input):
    if not isinstance(input,np.ndarray):
        raise TypeError(('Only numpy arrays are allowed to'
                         +'represent 3-component vectors'
                         +' not {}'.format(type(input))))
    if input.ndims != 2 or input.shape[1] != 3:
        raise ValueError(('Only [n x 3] numpy arrays are allowed'
                          +'as representation of 3-component vectors'
                          +' not {}'.format(input.shape)))
