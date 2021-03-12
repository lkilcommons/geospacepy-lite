# (C) 2020 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam M. Kilcommons
import numpy as np
import traceback, logging,textwrap
from functools import wraps
import datetime
import bisect

log = logging.getLogger('dmsp.'+__name__)
#Prefixing with dmsp. will make log calls propagate up to a master logger called dmsp
#this makes it so we don't have to set logging settings here, just in the master logger

#Custom Libraries
import sys
#from geospacepy import lmk_utils #my toolbox of convenience functions

#J2000 Epoch (remember that J2000 is defined as starting at NOON on 1-1-2000 not midnight)
dt_j2000 = datetime.datetime(2000,1,1,12,0,0)

def _update_vectorized_conversion_func_docstring(vectorized_conversion_func):
    #First line is short description in numpy style docstrings
    origdoc = vectorized_conversion_func.__doc__
    short_description = origdoc.splitlines()[0]
    newdoc = origdoc.replace(short_description,short_description+' (vectorized)')
    vectorized_conversion_func.__doc__ = newdoc
    
def _vectorized_datetime2(conversion_func):
    """Generates a vectorized version of a conversion function which converts
    from datetime to something else"""
    @wraps(conversion_func)
    def vectorized_conversion_func(datetimes,*args):

        if isinstance(datetimes,datetime.datetime):
            raise ValueError(('This functions is for lists/arrays of datetimes.'
                              +'Use {}'.format(conversion_func.__name__)
                              +'for single (scalar) datetime inputs'))
        elif isinstance(datetimes,np.ndarray):
            #Numpy array of datetimes passed
            dts = datetimes.flatten().tolist()
            converted_shape = datetimes.shape
            converted_size = datetimes.size
        else:
            #Make sure it's indexable
            try:
                dt = datetimes[0]
            except:
                raise ValueError(('Datetime collection {}'.format(datetimes)
                                  +' is not indexable. Try an array or list.'))
            #Assume things will work out and plan to return a flat array
            dts = datetimes
            converted_shape = (len(datetimes),)
            converted_size = len(datetimes)

        #Handle the arguments
        #Strings are indexable but not a collection, so forbid using them
        #as arguments to a time conversion function
        if any([isinstance(arg,str) for arg in args]):
            raise ValueError(('An argument of type str used in a time conversion'
                              +'function. Since strings are indexable and'
                              +'not easily distinguishable from "real"'
                              +'collections, this is not an allowed input type'))

        cleaned_args = []
        arg_is_indexable = []
        for iarg,arg in enumerate(args):
            if isinstance(arg,np.ndarray):
                #Passed argument is a np array
                if len(arg)==1 or arg.size==converted_size:
                    cleaned_arg = arg.flatten().tolist()
                    if len(arg)==1:
                        cleaned_args.append(cleaned_arg[0])
                        arg_is_indexable.append(False)
                    else:
                        cleaned_args.append(cleaned_arg)
                        arg_is_indexable.append(True)
                else:
                    raise ValueError(('Extra argument {} '.format(iarg)
                                      +'unexpected shape '
                                      +' ({}!={})'.format(arg.shape,
                                                        converted_shape)))
            elif isinstance(arg,float) or isinstance(arg,int):
                #A scalar value was passed
                cleaned_args.append(arg)
                arg_is_indexable.append(False)
            else:
                try:
                    a = arg[0]
                    arglen = len(arg)
                except:
                    raise ValueError(('Argument {} is not indexable'.format(arg)
                                      +' or does not have a length'
                                      +' if it is meant to be a collection'
                                      +' of inputs to a time conversion function'
                                      +' try using an array or list'))
                cleaned_args.append(arg)
                arg_is_indexable.append(True)

        #Finally do the conversion
        converted_times = np.empty(converted_shape)
        for i_t in range(len(dts)):
            conversion_args = [dts[i_t]]
            for arg,indexable in zip(cleaned_args,arg_is_indexable):
                if indexable:
                    conversion_args.append(arg[i_t])
                else:
                    conversion_args.append(arg)

            #This is required to make outputing column arrays work
            ind_tup = np.unravel_index(i_t,converted_shape)
            converted_times[ind_tup] = conversion_func(*conversion_args)

        return converted_times

    vectorized_name = conversion_func.__name__.replace('datetime2','datetimearr2')
    vectorized_conversion_func.__name__ = vectorized_name

    _update_vectorized_conversion_func_docstring(vectorized_conversion_func)
    
    return vectorized_conversion_func

def _vectorized_2datetime(conversion_func):
    """Generates a function which converts a time into an
    np array of datetimes. If the inputs are np arrays, the output array
    will mirror the inputs' shape. If inputs aren't arrays, then the output
    array will be flat"""
    @wraps(conversion_func)
    def vectorized_conversion_func(numerical_times,*args):
        #Strings are indexable but not a collection, so forbid using them
        #as arguments to a time conversion function
        if any([isinstance(arg,str) for arg in args]):
            raise ValueError(('An argument of type str used in a time conversion'
                              +'function.'))

        if isinstance(numerical_times,float) or isinstance(numerical_times,int):
            raise ValueError(('This function is for lists/arrays'
                              +'Use {}'.format(conversion_func.__name__)
                              +'for single (scalar) inputs'))

        elif isinstance(numerical_times,np.ndarray):
            #Numpy array of datetimes passed
            ts = numerical_times.flatten().tolist()
            converted_shape = numerical_times.shape
            converted_size = numerical_times.size
        else:
            #Make sure it's indexable and has a length
            try:
                t = numerical_times[0]
                tlen = len(numerical_times)
            except:
                raise ValueError(('Time input {}'.format(numerical_times)
                                  +' is not indexable or does not'
                                  +' have a length. Try an array or list.'))
            #Assume things will work out and plan to return a flat array
            ts = numerical_times
            converted_shape = (len(ts),)
            converted_size = len(ts)

        #Handle the arguments
        #Strings are indexable but not a collection, so forbid using them
        #as arguments to a time conversion function
        if any([isinstance(arg,str) for arg in args]):
            raise ValueError(('An argument of type str used in a time conversion'
                              +'function. Since strings are indexable and'
                              +'not easily distinguishable from "real"'
                              +'collections, this is not an allowed input type'))

        #Check the rest of the arguments
        cleaned_args = []
        arg_is_indexable = []
        for iarg,arg in enumerate(args):
            if isinstance(arg,np.ndarray):
                #Passed argument is a np array
                if len(arg)==1 or arg.size==converted_size:
                    cleaned_arg = arg.flatten().tolist()
                    if len(arg)==1:
                        cleaned_args.append(cleaned_arg[0])
                        arg_is_indexable.append(False)
                    else:
                        cleaned_args.append(cleaned_arg)
                        arg_is_indexable.append(True)
                else:
                    raise ValueError(('Extra argument {} '.format(iarg)
                                      +'unexpected shape '
                                      +' ({}!={})'.format(arg.shape,
                                                        converted_shape)))
            elif isinstance(arg,float) or isinstance(arg,int):
                #A scalar value was passed
                cleaned_args.append(arg)
                arg_is_indexable.append(False)
            else:
                try:
                    a = arg[0]
                    arglen = len(arg)
                except:
                    raise ValueError(('Argument {} is not indexable'.format(arg)
                                      +' or does not have a length'
                                      +' if it is meant to be a collection'
                                      +' of inputs to a time conversion function'
                                      +' try using an array or list'))
                cleaned_args.append(arg)
                arg_is_indexable.append(True)

        #Initialize the output array
        converted_dts = np.empty(converted_shape,dtype='object')
        print(converted_shape)
        for i_t in range(len(ts)):
            conversion_args = [ts[i_t]]
            for arg,indexable in zip(cleaned_args,arg_is_indexable):
                if indexable:
                    conversion_args.append(arg[i_t])
                else:
                    conversion_args.append(arg)
            #This is required to make outputing column arrays work
            ind_tup = np.unravel_index(i_t,converted_shape)
            converted_dts[ind_tup] = conversion_func(*conversion_args)

        return converted_dts

    vectorized_name = conversion_func.__name__.replace('2datetime','arr2datetime')
    vectorized_conversion_func.__name__ = vectorized_name

    _update_vectorized_conversion_func_docstring(vectorized_conversion_func)
   
    return vectorized_conversion_func

def datetime2jd(dt):
    """Converts between Python datetime and Julian Date
    (days since 12:00 PM on January 1, 4713 B.C.)
    Implementation is valid from 1900-2199

    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    jd : float
    """
    if dt.year < 1900:
        raise ValueError('Year must be 4 digit year')
    t1 = 367.*dt.year
    t2 = np.floor(7.*(dt.year+np.floor((dt.month+9.)/12.))/4.)
    t3 = np.floor(275.*dt.month/9.)
    t4 = dt.day + 1721013.5
    t5 = (((dt.microsecond/1e6+dt.second)/60.+dt.minute)/60+dt.hour)/24
    jd = t1-t2+t3+t4+t5
    return jd

datetimearr2jd = _vectorized_datetime2(datetime2jd)

def jd2datetime(jd):
    """Converts between Julian Date (days since 12:00 PM on January 1, 4713 B.C.)
    and Python datetime.
    Implementation is valid for 1900-2199

    Parameters
    ----------
    jd : float

    Returns
    -------
    dt : datetime.datetime

    """
    T1900 = (jd-2415019.5)/365.25
    if T1900<0:
        raise ValueError('Cannot convert Julian dates before 12:00 Jan 1 1900')
    year = 1900+np.floor(T1900)
    leapyrs = np.floor((year-1900-1)*.25)
    days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    if days < 1.0:
        year-=1
        leapyrs = np.floor((year-1900-1)*.25)
        days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    #doy = np.floor(days)
    return datetime.datetime(int(year),1,1,0,0,0)+datetime.timedelta(days=(days-1))

jdarr2datetime = _vectorized_2datetime(jd2datetime)

def datetime2j2000(dt):
    """Datetime to Julian date relative to j2000 Epoch (Noon on Jan 1, 2000)
    
    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    j2000 : float
    
    """
    return datetime2jd(dt)-datetime2jd(dt_j2000)

def j20002datetime(j2000):
    """Julian date relative to j2000 Epoch (Noon on Jan 1, 2000) to datetime
    
    Parameters
    ----------
    j2000 : float

    Returns
    -------
    dt : datetime.datetime
    """
    return jd2datetime(j2000 + datetime2jd(dt_j2000))

def datetime2doy(dt):
    """Python datetime to (decimal) day of year
    
    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    doy : float
        
    """
    return dt.timetuple().tm_yday + dt.hour/24. + dt.minute/24./60. + dt.second/86400. + dt.microsecond/86400./1e6

datetimearr2doy = _vectorized_datetime2(datetime2doy)

def doy2datetime(doy,year):
    """Day of year to python datetime
    
    Parameters
    ----------
    doy : float
    year : float or int

    Returns
    -------
    dt : datetime.datetime
    """
    return datetime.datetime(int(year),1,1,0,0,0)+datetime.timedelta(days=doy-1.) #Returns floating point day of year

doyarr2datetime = _vectorized_2datetime(doy2datetime)

def datetime2datenum(dt):
    """Python datetime to matlab epoch (datenum)

    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    mlep : float
    """
    return dt.toordinal() + 366. + dt.hour/24. + dt.minute/24./60. + dt.second/86400. + dt.microsecond/86400./1e6

datetimearr2datenum = _vectorized_datetime2(datetime2datenum)

def datenum2datetime(mlep):
    """Matlab epoch (datenum) to python datetime
    
    Parameters
    ----------
    mlep : float

    Returns
    -------
    dt : datetime.datetime
    """
    return datetime.datetime.fromordinal(int(np.floor(mlep))) + datetime.timedelta(days=np.mod(mlep,1)) - datetime.timedelta(days = 366)

datenumarr2datetime = _vectorized_2datetime(datenum2datetime)

def datetime2sod(dt):
    """Python datetime to (decimal) second of day

    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    sod : float
    """
    return (dt-datetime.datetime.combine(dt.date(),datetime.time(0))).total_seconds()

datetimearr2sod = _vectorized_datetime2(datetime2sod) #Not checked

def sod2datetime(sod,year,month,day):
    """Second of day to datetime
    
    Parameters
    ----------
    sod : float
    year : float or int
    month : float or int
    day : float or int

    Returns
    -------
    sod : float
    """
    return datetime.datetime(year,month,day)+datetime.timedelta(seconds=sod)

sodarr2datetime = _vectorized_2datetime(sod2datetime)

def datetime2soy(dt):
    """Datetime to second of year
    
    Parameters
    ----------
    dt : datetime.datetime

    Returns
    -------
    soy : float
    """
    return (dt-datetime.datetime(dt.year,1,1)).total_seconds()

def soy2datetime(soy,year):
    """Second of year to datetime
    
    Parameters
    ----------
    soy : float
    year : float or int

    Returns
    -------
    dt : datetime.datetime
    """
    return datetime.datetime(year,1,1)+datetime.timedelta(seconds=np.floor(soy))

# def datetimearr2jd(datetimearr):
#     """
#     Converts a n x 1 or 1 x n or (n,) array of python datetimes
#     to an n x 1 julian day array
#     """
#     #Input sanitize
#     if isinstance(datetimearr,np.ndarray):
#         datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

#     jd = np.empty((len(datetimearr),1))
#     for k,dt in enumerate(datetimearr):
#         jd[k,0] = datetime2jd(dt)
#     return jd

# def jdarr2datetime(jdarr):
#     """
#     Converts a n x 1 or 1 x n or (n,) array of julian days
#     to an n x 1 python datetime array
#     """
#     #Input sanitize
#     if isinstance(jdarr,np.ndarray):
#         jdarr = jdarr.flatten().tolist() #Easier to deal with list, no indexing ambiguity

#     dts = np.empty((len(jdarr),1),dtype='object')
#     for k in range(len(jdarr)):
#         dts[k,0] = jd2datetime(jdarr[k])

#     return dt

# def datetimearr2j2000(datetimearr):
#     """
#     Converts a n x 1 or 1 x n or (n,) array of python datetimes
#     to an n x 1 j2000 epoch julian date array
#     """
#     #Input sanitize
#     if isinstance(datetimearr,np.ndarray):
#         datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

#     j2000 = np.empty((len(datetimearr),1))
#     for k,dt in enumerate(datetimearr):
#         j2000[k,0] = datetime2j2000(dt)
#     return j2000

# def datetimearr2soy(datetimearr):
#     """
#     Converts a n x 1 array of python datetimes
#     to an n x 1 day of year array
#     """
#     #Input sanitize
#     if isinstance(datetimearr,np.ndarray):
#         datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

#     soy = np.empty((len(datetimearr),1))
#     for k,dt in enumerate(datetimearr):
#         soy[k,0] = datetime2soy(dt)
#     return soy

# def soyarr2datetime(soyarr,year):
#     """
#     Converts an n x 1 array of days of year
#     and either a single year value or a n x 1 array of years
#     to a np array of datetimes
#     """
#     #Input sanitize - make sure everything is a list
#     soyarr = lmk_utils.arraycheck(soyarr,returnList=True)
#     year = lmk_utils.arraycheck(year,returnList=True)

#     dt = np.empty((len(soyarr),1),dtype='object')
#     for k in range(len(soyarr)):
#         if len(year)>1:
#             dt[k,0] = soy2datetime(soyarr[k],np.floor(year[k,0]))
#         else:
#             dt[k,0] = soy2datetime(soyarr[k],np.floor(year[0][0]))

#     return dt

# def datetimearr2doy(datetimearr):
#     """
#     Converts a n x 1 array of python datetimes
#     to an n x 1 day of year array
#     """
#     #Input sanitize
#     if isinstance(datetimearr,np.ndarray):
#         datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

#     doy = np.empty((len(datetimearr),1))
#     for k,dt in enumerate(datetimearr):
#         doy[k,0] = datetime2doy(dt)
#     return doy

# def doyarr2datetime(doyarr,year):
#     """
#     Converts an n x 1 array of days of year
#     and either a single year value or a n x 1 array of years
#     to a np array of datetimes
#     """
#     #Input sanitize - make sure everything is a list
#     doyarr = lmk_utils.arraycheck(doyarr,returnList=True)
#     year = lmk_utils.arraycheck(year,returnList=True)

#     dt = np.empty((len(doyarr),1),dtype='object')
#     for k in range(len(doyarr)):
#         if len(year)>1:
#             dt[k,0] = doy2datetime(doyarr[k],np.floor(year[k]))
#         else:
#             dt[k,0] = doy2datetime(doyarr[k],np.floor(year[0][0]))

#     return dt

# def datetimearr2sod(datetimearr):
#     """
#     Converts a n x 1 array of python datetimes
#     to an n x 1 second of day array
#     """
#     #Input sanitize
#     if isinstance(datetimearr,np.ndarray):
#         datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

#     sod = np.empty((len(datetimearr),1))
#     for k,dt in enumerate(datetimearr):
#         sod[k,0] = datetime2sod(dt)
#     return sod

# def datearr2sod(datearr):
#     """
#     Converts a n x 6 array with columns year,month,day,hour,minute,second
#     to an n x 1 second of day array
#     This and the next method show two ways of dealing with casting
#     the input to int. Which is faster?
#     """
#     #Input sanitize
#     if not isinstance(datearr,np.ndarray):
#         raise RuntimeError('Input must be np.ndarray n x 6')

#     sod = np.empty((len(datearr[:,0]),1))
#     for k in range(len(datearr[:,0])):
#         sod[k,0] = ymdhms2sod(np.floor(datearr[k,0]),np.floor(datearr[k,1]),np.floor(datearr[k,2]),np.floor(datearr[k,3]),np.floor(datearr[k,4]),np.floor(datearr[k,5]))
#     return sod

# def sodarr2datetime(sodarr,year,month,day):
#     """
#     Converts a nx1 array of seconds of the day
#     to a datetime array. Year, month, day are
#     scalar integers. Returns a np array of
#     datetime objects.
#     """

#     #Input sanitize
#     if isinstance(sodarr,np.ndarray):
#         sodarr = sodarr.flatten().tolist() #Easier to deal with list, no indexing ambiguity

#     dt = np.empty((len(sodarr),1),dtype='object')
#     for k in range(len(sodarr)):
#         dt[k,0] = sod2datetime(sodarr[k],year,month,day)
#     return dt

# def datearr2doy(datearr):
#     """
#     Converts a n x 6 array with columns year,month,day,hour,minute,second
#     to an n x 1 decimal day of year array
#     """

#     #Input sanitize
#     if not isinstance(datearr,np.ndarray):
#         raise RuntimeError('Input must be np.ndarray n x 6')

#     datearr = datearr.astype('int')
#     doy = np.array((len(datearr[:,0]),1))
#     for k in range(len(datearr[:,0])):
#         doy[k,0] = ymdhms2doy(datearr[k,0],datearr[k,1],datearr[k,2],datearr[k,3],datearr[k,4],datearr[k,5])
#     return doy

# def datearr2datetime(datearr,asArray=True):
#     """
#     Converts a n x 6 array with columns year,month,day,hour,minute,second
#     to an n x 1 list, or if asArray==True, as a np array
#     This is sort of slow, but it gets the job done
#     """
#     datetimes = []
#     datearr = datearr.astype('int')
#     for k in range(len(datearr[:,0])):
#         datetimes.append(datetime.datetime(datearr[k,0],datearr[k,1],datearr[k,2],datearr[k,3],datearr[k,4],datearr[k,5]))
#     if asArray:
#         datetimes = np.array(datetimes,ndmin=2)

#     return datetimes

# def datenumarr2datetime(datenumarr,asArray=True):
#     """
#     Converts a n x 1 array of matlab datenumbers to
#     to an n x 1 list, or if asArray==True, as a np array
#     This is sort of slow, but it gets the job done
#     """

#     datetimes = [datenum2datetime(datenumarr[k]) for k in range(len(datenumarr))]
#     if asArray:
#         datetimes = np.array(datetimes)
#     return datetimes

def fastMatchTimes(primary_dt,dt,tol_us=4e5,fail_on_duplicates=True,allow_duplicates=False):
    """Finds a matching timestamp in primary_dt
    within tolerance tol_us (given in microseconds)
    for every value in dt.
    
    Parameters
    ----------

        dt - np.array(dtype=object), size=n
            The array of timestamps (Python datetime) 
            for which you want to find matches 
        
        primary_dt - np.array(dtype=object), size=m
            The array of timestamps (Python datetimes) which will be searched

    Returns
    -------
    
        inds - np.ndarray, size=n
            Array of len(dt) of indices into
            primary_dt or NaN if no match was found
            for a partiuclar dt value.
    
    """

    import bisect

    #Formatting String
    fmstr = '%Y-%m-%d %H:%M:%S.%f'

    reftime = datetime.datetime.utcfromtimestamp(0)

    # Input sanitize
    if isinstance(dt,np.ndarray):
        dt = dt.flatten().tolist()
    if isinstance(primary_dt,np.ndarray):
        primary_dt = primary_dt.flatten().tolist()

    #Convert to unix epoch (seconds since hour 0 Jan 1 1970)
    log.debug("Computing UNIX epoch milisecond integer values for all datetimes to improve matching speed")
    dt = np.array([(d-reftime).total_seconds() for d in dt])
    primary_dt = np.array([(d-reftime).total_seconds() for d in primary_dt])

    # Throw an error if there are NaN values
    #if any(np.isnan(dt)) or any(np.isnan(primary_dt)):
    #   raise RuntimeError('NaN datetimes found! Sanitize your inputs before calling matchTimes!')
    #Nope...apparently isnan doesn't work on datetime arrays

    #Preallocate the results array and fill it with nan
    inds = np.zeros((len(dt),1))
    inds.fill(np.nan)
    #inds = inds.flatten().tolist()

    starttime = datetime.datetime.now()

    ep2dt = lambda ep: reftime + datetime.timedelta(seconds=ep)

    log.debug("Beginning matching")
    n_unmatched = 0
    for i in range(len(dt)):
        if np.mod(i,5000)==0 and i>0:
            rate = float(i)/(datetime.datetime.now()-starttime).total_seconds()
            timetoend = (len(dt)-i)/rate/60.
            log.info("%d/%d matched, %d unmatched: ~%.3f iterations/sec. About %.1f minutes remain." % (i,len(dt),n_unmatched,rate,timetoend))
            n_unmatched = 0
        try:
            #Find the place where dt[i,0] should be inserted into primary_dt to maintain
            #sorted order. bisect and bisect_left should always return the same index
            #UNLESS the exact value of dt[i,0] is present
            ind = np.searchsorted(primary_dt,dt[i],side='right')

            if ind < len(primary_dt) and dt[i] == primary_dt[ind]: #This means that the exact value of dt[i] is present at ind
                inds[i] = ind

            elif ind > 0 and ind < len(primary_dt): #this means that primary_dt[ind] > dt[i] and primary_dt[ind-1] < dt[i]
                ind_l=ind-1

                #Compute deltas and compare to tolerence
                #delta_before = dt[i,0]-primary_dt[ind-1,0]
                #delta_after = primary_dt[ind,0]-dt[i,0]

                #Compute before and after time difference in microseconds
                delta_before_us = abs((dt[i]-primary_dt[ind-1]))*1.0e6
                delta_after_us = abs((primary_dt[ind]-dt[i]))*1.0e6

                #Determine if either of the values we'd like to assign is already matched
                before_duplicate = True if ind-1 in inds else False
                after_duplicate = True if ind in inds else False

                #Handle case of ambiguous match
                #Default behaviour is to choose best match if Ambigious

                #Determine if the time difference to the index before/after is within tolerence
                if delta_after_us < tol_us and delta_before_us < tol_us:

                    #Determine which value provides the smallest difference in timestamp
                    #best_match = 'before' if delta_before_us<=delta_after_us else 'after'

                    #Believe it or not, I've had the case happen where delta_before exactly equals delta_after. In
                    #this case, we assume that we want to use the earlier timestamp because of 'the arrow of time'
                    #i.e. using information from the past is less uncertain that from the future, all other things
                    #being equal

                    log.warn('Ambiguous match between two timestamps.\n%s within tolerance %f s of both\n %s (%f s) and\n %s (%f s)' % (ep2dt(dt[i]).strftime(fmstr),tol_us/1.0e6,\
                            ep2dt(primary_dt[ind]).strftime(fmstr),delta_after_us/1.0e6,
                            ep2dt(primary_dt[ind-1]).strftime(fmstr),delta_before_us/1.0e6) )
                    if delta_before_us<=delta_after_us and (not before_duplicate or allow_duplicates):
                        inds[i] = ind-1
                    elif delta_before_us>delta_after_us and (not after_duplicate or allow_duplicates):
                        inds[i] = ind
                else:
                    #No need to use best_match here, we only have one value within tolerance
                    if delta_before_us < tol_us and (not before_duplicate or allow_duplicates):
                        inds[i] = ind-1
                    elif delta_after_us < tol_us and (not after_duplicate or allow_duplicates):
                        inds[i] = ind
                    else:
                        n_unmatched += 1
                        #log.warn('No match found for timestamp %s' % (ep2dt(dt[i]).strftime(fmstr)))

            elif ind == 0 or ind == len(primary_dt):
                #Edge cases, use bisect left result if within tolerence
                if ind==len(primary_dt):
                    ind-=1
                delta = dt[i] - primary_dt[ind]
                if abs(delta*1.0e6) < tol_us and (ind not in inds or allow_duplicates):
                    inds[i] = ind
                else:
                    n_unmatched += 1
                    #log.warn('No match found for timestamp %s [NOTE: Edge case, closest match was beginning or end of primary_dt]' % \
                    #   (str(dt[i])))

            else:
                log.error('Unhandled bisect case...this should not happen')
                raise RuntimeError('Something wierd is going on, unhandled bisect case')

        except:
            traceback.print_exc()
            print("Error in matchTimes:\nIndex in dt: %d\nBisect Right Index in primary_dt: %d\nBisect Left Index in primary_dt: %d" % (i,ind,ind_l))
            print("dt[i] %s\nprimary_dt[ind] %s\nprimary_dt[ind_l] %s" % (str(dt[i]),str(primary_dt[ind]),str(primary_dt[ind_l])))
            #pdb.set_trace()

    #print "Completed timestamp alignment, %d points remained unmatched out of %d\n" %(unmatched_counter,len(dt))

    #Convert list back to array
    #inds = np.array(inds)
    inds = inds.astype(int)
    #inds = lmk_utils.asColumn(inds)

    #Identify the first occurance of each unique value in inds in unqinds
    unqvals,unqinds = np.unique(inds,return_index=True)
    #Identifiy all indices in inds which are nan
    naninds = np.nonzero(np.isnan(inds))

    nunq = len(unqvals)
    log.info("%d/%d values to match were matched uniquely." % (nunq,len(inds)))
    log.info("%d/%d values are nan (unmatched)." % (len(naninds),len(inds)))

    return inds
