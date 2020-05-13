# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import numpy as np
import inspect
import functools

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
        l = number.size
        if l != 1:
            return False
        else:
            return True
    except AttributeError:
        if isinstance(number,int) or isinstance(number,float):
            return True
        else:
            return False

def BroadcastLenOneInputsToMatchArrayInputs(wrapped_func):
    """Decorator which broadcasts any single valued input 
    (numerical primative or array of length 1) into a numpy array
    which matches the size of the other inputs before passing
    the inputs to the wrapped function"""
    @functools.wraps(wrapped_func)
    def wrapper(*args,**kwargs):
        args_shapes = []
        unique_shapes = []
        for arg in args:
            if is_number_or_len_one_array(arg):
                arg_shape = (1,)
            else:
                arg_shape = arg.shape
            args_shapes.append(arg_shape)
            if arg_shape not in unique_shapes:
                unique_shapes.append(arg_shape)

        unique_nonsingle_shapes = [shape for shape in unique_shapes if shape!=(1,)]

        if len(unique_nonsingle_shapes)>1:
            raise ValueError(('Unable to broadcast inputs '
                              +'of shapes {} '.format(args_shapes)
                              +'into a single common shape '
                              +'wrapped function is {}'.format(wrapped_func)))
        if len(unique_nonsingle_shapes)==0:
            common_shape = (1,)
        else:
            common_shape = unique_nonsingle_shapes[0]
            
        broadcasted_args = []
        for arg,shape in zip(args,args_shapes):
            if shape == common_shape:
                broadcasted_args.append(arg)
            else:
                broadcasted_args.append(np.ones(common_shape)*arg)
        
        return wrapped_func(*broadcasted_args,**kwargs)
    return wrapper

def _check_follows_3_component_vector_convention(inpt):
    """Determine if an given object is a numpy array with n rows and
    3 columns (geospacepy-lite convention for a 3 component vector)

        PARAMETERS
        ----------

            inpt - object
                The object (of any type) to check

        RAISES
        ------

            TypeError
                If it's not an array
            ValueError
                If the shape is not (n,3)

    """
    if not isinstance(inpt,np.ndarray):
        raise TypeError(('Only numpy arrays are allowed to'
                         +' represent 3-component vectors'
                         +' not {}'.format(type(inpt))))
    if inpt.ndim != 2 or inpt.shape[1] != 3:
        raise ValueError(('Only [n x 3] numpy arrays are allowed'
                          +' as representation of 3-component vectors'
                          +' not {}'.format(inpt.shape)))

class CheckInputsAreThreeComponentVectors(object):
    """Decorator which type and shape
    checks the wrapped function's inputs to ensure they follow the 
    shape=(n,3) convention for 3 component vectors"""  
    def __init__(self,vector_args_names):
        """

        PARAMETERS
        ----------

        vector_args_names : str or list
            Names of the wrapped function's
            inputs which will be type-and-shape checked
        """
        if not isinstance(vector_args_names,list):
            self.vector_args_names = [vector_args_names]
        else:
            self.vector_args_names = vector_args_names

    @staticmethod
    def _determine_argument_names(func):
        """Find the names of the named, non-keyword arguments of function
        func (ignores input=default keyword inputs and *args, **kwargs inputs)
        """
        sig = inspect.signature(func)
        func_args = []
        for pname,param in sig.parameters.items():
            #Check if argument could be from *args or **kwargs
            #(it's kind would be VAR_POSITIONAL or VAR_KEYWORD)
            is_not_star = param.kind==inspect.Parameter.POSITIONAL_OR_KEYWORD
            is_not_kwarg = param.default is param.empty
            if is_not_star and is_not_kwarg:
                func_args.append(pname)
        return func_args

    def _check_passed_argnames_are_valid_func_argnames(self,func):
        func_argnames = self._determine_argument_names(func)
        for vec_argname in self.vector_args_names:
            if vec_argname not in func_argnames:
                raise ValueError(('Function {}'.format(func.__name__)
                                +'has no input {}'.format(vec_argname)))

    def __call__(self,func):
        self._check_passed_argnames_are_valid_func_argnames(func) 
        func_argnames = self._determine_argument_names(func)
        argnums_to_check = []
        for argname in self.vector_args_names:
            argnums_to_check.append(func_argnames.index(argname))
        
        @functools.wraps(func) 
        def wrapper(*args,**kwargs):
            for iarg in argnums_to_check:
                _check_follows_3_component_vector_convention(args[iarg])
            return func(*args,**kwargs)
        
        return wrapper

