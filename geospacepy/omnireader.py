# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons

# This is a dummy module file to redirect users to the new omnireader package
class omni_interval(object):

    def __init__(self,*args,**kwargs):
        raise RuntimeError(('Omnireader is no longer part of Geospacepy-lite'
                            +'it has been moved to its own package at'
                            +' github.com/lkilcommons/nasaomnireader'))