import pytest
import numpy as np
import numpy.testing as nptest
from geospacepy.rotations import rot
from scipy.spatial.transform import Rotation

@pytest.mark.parametrize('n_vecs,n_angles',[(1,1),(3,1),(3,3)])
def test_90_degree_rotation_of_xhat_about_z_produces_yhat(n_vecs,n_angles):
    vecs = np.array([[1.,0.,0.] for i_vec in range(n_vecs)])
    expected_rotated_vecs = np.array([[0.,1.,0.] for i_vec in range(n_vecs)])
    angle = np.radians(-90.)
    angles = angle if n_angles == 1 else np.ones((n_angles,))*angle
    axis=2
    rotated_vecs = rot(angles,axis,vecs)
    nptest.assert_allclose(rotated_vecs,expected_rotated_vecs,atol=1e-8,rtol=0)

