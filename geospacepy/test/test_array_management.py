import pytest
import numpy as np
import numpy.testing as nptest
from geospacepy.array_management import (is_number_or_len_one_array,
                                        BroadcastLenOneInputsToMatchArrayInputs,
                                        CheckInputsAreThreeComponentVectors)

def test_float_is_number_or_len_one_array():
    assert is_number_or_len_one_array(1.)

@pytest.mark.parametrize('arr',(np.array([1.]),np.array([[1.]])))
def test_arr_is_number_or_len_one_array(arr):
    assert is_number_or_len_one_array(arr)

@pytest.mark.parametrize('arr',(np.array([1.,2.]),np.array([[1.],[2.]])))
def test_arr_is_NOT_number_or_len_one_array(arr):
    assert not is_number_or_len_one_array(arr)

def _cart2sph(x,y,z):
    """A function that should operate equally well on floats and arrays,
    and involves trignometry...a good test function for the types of
    functions in geospacepy-lite"""
    r = x**2+y**2+z**2
    th = np.arctan2(y,x)
    ph = np.arctan2(x**2+y**2,z)
    return r,th,ph

def _cart2sph_example():
    x,y,z = 0.,1.,0.    
    r,th,ph = 1.,np.pi/2.,np.pi/2.
    return x,y,z,r,th,ph 

def test_cart2sph_example():
    """Sanity check the example"""
    x,y,z,r_ref,th_ref,ph_ref = _cart2sph_example()
    r,th,ph = _cart2sph(x,y,z)
    tol = .0001
    assert(np.abs(r-r_ref)<tol)
    assert(np.abs(th-th_ref)<tol)
    assert(np.abs(ph-ph_ref)<tol)

def _bcast_to_shape(value,shape):
    if shape == 1:
        return value
    else:
        return np.ones(shape)*value
    
@pytest.mark.parametrize('xshp,yshp,zshp,expctshp',[(1,1,1,1),
                                            ((3,),1,1,(3,)),
                                            (1,(1,3),(1,3),(1,3))])
def test_broadcast_len_one_inputs_to_match_array_inputs_on_cart2sph(xshp,yshp,zshp,
                                                                    expctshp):
    """Test the BroadcastLenOneInputsToMatchArrayInputs decorator
    by calling a decorated function (which converts cartesian to spherical coords)
    """
    wrapped_cart2sph = BroadcastLenOneInputsToMatchArrayInputs(_cart2sph)

    x,y,z,r,th,ph = _cart2sph_example()
    xarr = _bcast_to_shape(x,xshp)
    yarr = _bcast_to_shape(y,yshp)
    zarr = _bcast_to_shape(z,zshp)

    rarr,tharr,pharr = wrapped_cart2sph(xarr,yarr,zarr)

    tol = .0001
    if expctshp == 1:
        assert(np.abs(rarr-r)<tol)
        assert(np.abs(tharr-th)<tol)
        assert(np.abs(pharr-ph)<tol)
    else:
        assert(rarr.shape==expctshp)
        assert(tharr.shape==expctshp)
        assert(pharr.shape==expctshp)

def _vec_cart2sph(Rxyz):
    x=Rxyz[:,0].flatten()
    y=Rxyz[:,1].flatten()
    z=Rxyz[:,2].flatten()
    r,th,ph = _cart2sph(x,y,z)
    Rrthph = np.column_stack([r,th,ph])
    return Rrthph

def _vec_cart2sph_example():
    x,y,z,r,th,ph = _cart2sph_example()
    Rxyz = np.array([x,y,z]).reshape(-1,3)
    Rrthph = np.array([r,th,ph]).reshape(-1,3)
    return Rxyz,Rrthph

@CheckInputsAreThreeComponentVectors('Rxyz')
def _decorated_vec_cart2sph(Rxyz):
    return _vec_cart2sph(Rxyz)

@pytest.mark.parametrize('n_vectors',[1,3])
def test_check_inputs_are_three_component_vectors_on_decorated_vec_cart2sph(n_vectors):
    Rxyz,Rrthph = _vec_cart2sph_example()

    Rxyz_in = np.broadcast_to(Rxyz,(n_vectors,3))
    Rrthph_ref = np.broadcast_to(Rrthph,(n_vectors,3))

    Rrthph_out = _decorated_vec_cart2sph(Rxyz_in)

    nptest.assert_allclose(Rrthph_out,Rrthph_ref)

def test_check_inputs_are_three_component_vectors_raises_type_error_on_nonarray():
    with pytest.raises(TypeError):
        not_a_vector = 'reallyimavector'
        _decorated_vec_cart2sph(not_a_vector)

def test_check_inputs_are_three_component_vectors_raises_value_error_on_bad_shape():
    with pytest.raises(ValueError):
        Rxyz,Rrthph = _vec_cart2sph_example()
        _decorated_vec_cart2sph(Rxyz.T)
