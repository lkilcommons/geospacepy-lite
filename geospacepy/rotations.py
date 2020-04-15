import numpy as np

from geospacepy.array_management import is_number_or_len_one_array

def rotmat(angle,axis):
    """Return the 3x3 elementary rotation matrix (M),
    which when pre-multiplied (M*v) by a 3x1 column vector (v)
    returns the new components of v in a coordinate system 
    which has been rotated through an angle (angle) in *radians*
    about axis (axis), relative to the vector's native coorinate system
    Note: Positive angles represent rotations according the the right hand rule 
    (the directions of the fingers if the thumb is the axis)
    """ 

    if not is_number_or_len_one_array(angle):
        raise TypeError('Invalid angle value {}'.format(angle))

    if axis not in [0,1,2]:
        raise ValueError('Invalid value for axis, valid values are 0,1 and 2')

    c,s = np.cos(angle),np.sin(angle)
    if axis == 0:
        M = np.array([[1,   0,  0],
                      [0,   c,  s],
                      [0,-1*s,  c]])
    elif axis == 1:
        M = np.array([[c,   0,-1*s],
                      [0,   1,   0],
                      [s,   0,   c]])
    elif axis == 2:
        M = np.array([[   c,   s,  0],
                     [-1*s,   c,  0],
                     [   0,  0,  1]])
    return M

def rot(angles,axis,vecs):
    """Return representation of 3-component vectors vecs in a
    coordinate systems which has been rotated an angle (angle), about
    axis (axis)
    """

    if vecs.shape[1] != 3:
        raise ValueError(('Input vecs does have a dimension 1 of size 3'
                          'for multiple vectors, this function expects'
                          'a tall form input [n_vecs x 3]'))

    n_vecs = vecs.shape[0]
    if is_number_or_len_one_array(angles):
        angles_for_rot = angles*np.ones((n_vecs,))
    elif angles.size == n_vecs:
        angles_for_rot = angles.reshape((-1,1))
    else:
        raise ValueError('Cannot map {} angle to {} vectors'.format(angles.size,
                                                                    n_vecs))
    
    rotated_vecs = np.full_like(vecs, np.nan)
    for i_vec in range(n_vecs):
        M = rotmat(angles_for_rot[i_vec],axis)
        v = vecs[i_vec,:].T
        rotated_vecs[i_vec,:] = np.dot(M,v)

    return rotated_vecs

def rot1(angles,vecs):
    return rot(angles,0,vecs)

def rot2(angles,vecs):
    return rot(angles,1,vecs)

def rot3(angles,vecs):
    return rot(angles,2,vecs)

if __name__ == '__main__':

    from scipy.spatial.transform import Rotation

    vecs = np.array([[1,0,0],
                     [0,1,0],
                     [0,0,1]],dtype=float)

    angle = np.radians(45)

    rotated_vecs = rot(angle,2,vecs)
    print('Mine (1 angle) \n {}'.format(rotated_vecs))

    r = Rotation.from_euler('z',-1*angle)
    sp_rotated_vecs = r.apply(vecs)
    print('Scipy (1 angle) \n {}'.format(sp_rotated_vecs))

    angles = np.array([np.pi/2,np.pi,3*np.pi/2])
    rotated_vecs = rot(angles,2,vecs)
    print('Mine (3 angles) \n {}'.format(rotated_vecs))

    r = Rotation.from_euler('z',-1*angles)
    sp_rotated_vecs = r.apply(vecs)
    print('Scipy (3 angles) \n {}'.format(sp_rotated_vecs))
