"""
    special_datetime.py
    A library of miscellaneous useful vectorized datetime conversion functions.
    Also has timestamp matching routine
"""
import numpy, pdb, traceback, logging
import datetime
import bisect

log = logging.getLogger('dmsp.'+__name__)
#Prefixing with dmsp. will make log calls propagate up to a master logger called dmsp
#this makes it so we don't have to set logging settings here, just in the master logger

#Custom Libraries
import sys
from geospacepy import lmk_utils #my toolbox of convenience functions

#Date vector to second of day
def ymdhms2sod(y,mo,d,h,m,s):
    return float((datetime.datetime(y,mo,d,h,m,s)-\
    datetime.datetime.combine(datetime.date(y,mo,d),datetime.time(0))).total_seconds())

#Date vector to decimal day of year
def ymdhms2doy(y,mo,d,h,m,s):
    return float(datetime.datetime(y,mo,d,h,m,s).timetuple().tm_yday)+date2sod(y,mo,d,h,m,s)/86400.

#TODO Write Inverse: jd2datetime
def ymdhms2jd(year,mon,day,hr,mn,sc):
    #Takes UTC ymdhms time and returns julian date
    #FIXME: Add leap second support
    leapsecond = False
    if year < 1900:
        raise ValueError('Year must be 4 digit year')
    t1 = 367.*year

    t2 = int(7.*(year+int((mon+9.)/12.))/4.)
    t3 = int(275.*mon/9.)
    t4 = day + 1721013.5
    if not leapsecond:
        t5 = ((sc/60.+mn)/60+hr)/24
    else:
        t5 = ((sc/61.+mn)/60+hr)/24
    #print t1,t2,t3,t4,t5
    return t1-t2+t3+t4+t5

def jd2ymdhms(jd):
    dt = jd2datetime(jd)
    return dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second

def jd2datetime(jd):
    #Takes julian date and returns datetime.datetime in UTC

    #The inverse of above, from Vallado pp 208. (algorithm 22)
    T1900 = (jd-2415019.5)/365.25
    year = 1900+int(T1900)
    leapyrs = int((year-1900-1)*.25)
    days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    if days < 1.0:
        year-=1
        leapyrs = int((year-1900-1)*.25)
        days = (jd-2415019.5)-((year-1900)*(365.0) + leapyrs)
    #doy = int(days)
    return datetime.datetime(year,1,1,0,0,0)+datetime.timedelta(days=days-1)

#Julian Date (reasonably precise, but w/o leap seconds)
def datetime2jd(dt):
    return ymdhms2jd(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)

#J2000 Epoch (remember that J2000 is defined as starting at NOON on 1-1-2000 not 0 UT)
def datetime2j2000(dt):
    return ymdhms2jd(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)-ymdhms2jd(2000,1,1,12,0,0)

#Matlab epoch (datenum) to python datetime
def datenum2datetime(mlep):
    return datetime.datetime.fromordinal(int(mlep)) + datetime.timedelta(days=mlep%1) - datetime.timedelta(days = 366)

#Python datetime to day of year
def datetime2doy(dt):
    return dt.timetuple().tm_yday + dt.hour/24. + dt.minute/24./60. + dt.second/86400. + dt.microsecond/86400./1e6

#Python datetime to matlab epoch (datenum)
def datetime2datenum(dt):
    return dt.toordinal() + 366. + dt.hour/24. + dt.minute/24./60. + dt.second/86400. + dt.microsecond/86400./1e6

#Python datetime to second of day
def datetime2sod(dt):
    return (dt-datetime.datetime.combine(dt.date(),datetime.time(0))).total_seconds() #Microsecond precision

#Day of year to python datetime (full precision)
def doy2datetime(doy,year):
    return datetime.datetime(year,1,1,0,0,0)+datetime.timedelta(days=doy-1.) #Returns floating point day of year

#Second of day to datetime
def sod2datetime(sod,year,month,day):
    return datetime.datetime(year,month,day)+datetime.timedelta(seconds=sod)

#Second of year to datetime
def soy2datetime(soy,year):
    return datetime.datetime(year,1,1)+datetime.timedelta(seconds=int(soy))

#Datetime to second of year
def datetime2soy(dt):
    return (dt-datetime.datetime(dt.year,1,1)).total_seconds()


def soyarr2datetime(soyarr,year):
    """
    Converts an n x 1 array of days of year
    and either a single year value or a n x 1 array of years
    to a numpy array of datetimes
    """
    #Input sanitize - make sure everything is a list
    soyarr = lmk_utils.arraycheck(soyarr,returnList=True)
    year = lmk_utils.arraycheck(year,returnList=True)

    dt = numpy.empty((len(soyarr),1),dtype='object')
    for k in range(len(soyarr)):
        if len(year)>1:
            dt[k,0] = soy2datetime(soyarr[k],int(year[k,0]))
        else:
            dt[k,0] = soy2datetime(soyarr[k],int(year[0][0]))

    return dt

def datetimearr2soy(datetimearr):
    """
    Converts a n x 1 array of python datetimes
    to an n x 1 day of year array
    """
    #Input sanitize
    if isinstance(datetimearr,numpy.ndarray):
        datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

    soy = numpy.empty((len(datetimearr),1))
    for k,dt in enumerate(datetimearr):
        soy[k,0] = datetime2soy(dt)
    return soy

def doyarr2datetime(doyarr,year):
    """
    Converts an n x 1 array of days of year
    and either a single year value or a n x 1 array of years
    to a numpy array of datetimes
    """
    #Input sanitize - make sure everything is a list
    doyarr = lmk_utils.arraycheck(doyarr,returnList=True)
    year = lmk_utils.arraycheck(year,returnList=True)

    dt = numpy.empty((len(doyarr),1),dtype='object')
    for k in range(len(doyarr)):
        if len(year)>1:
            dt[k,0] = doy2datetime(doyarr[k],int(year[k]))
        else:
            dt[k,0] = doy2datetime(doyarr[k],int(year[0][0]))

    return dt

def datetimearr2doy(datetimearr):
    """
    Converts a n x 1 array of python datetimes
    to an n x 1 day of year array
    """
    #Input sanitize
    if isinstance(datetimearr,numpy.ndarray):
        datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

    doy = numpy.empty((len(datetimearr),1))
    for k,dt in enumerate(datetimearr):
        doy[k,0] = datetime2doy(dt)
    return doy

def datetimearr2jd(datetimearr):
    """
    Converts a n x 1 or 1 x n or (n,) array of python datetimes
    to an n x 1 julian day array
    """
    #Input sanitize
    if isinstance(datetimearr,numpy.ndarray):
        datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

    jd = numpy.empty((len(datetimearr),1))
    for k,dt in enumerate(datetimearr):
        jd[k,0] = datetime2jd(dt)
    return jd

def jdarr2datetime(jdarr):
    """
    Converts a n x 1 or 1 x n or (n,) array of julian days
    to an n x 1 python datetime array
    """
    #Input sanitize
    if isinstance(jdarr,numpy.ndarray):
        jdarr = jdarr.flatten().tolist() #Easier to deal with list, no indexing ambiguity

    dts = numpy.empty((len(jdarr),1),dtype='object')
    for k in range(len(jdarr)):
        dts[k,0] = jd2datetime(jdarr[k])

    return dts

def datetimearr2j2000(datetimearr):
    """
    Converts a n x 1 or 1 x n or (n,) array of python datetimes
    to an n x 1 j2000 epoch julian date array
    """
    #Input sanitize
    if isinstance(datetimearr,numpy.ndarray):
        datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

    j2000 = numpy.empty((len(datetimearr),1))
    for k,dt in enumerate(datetimearr):
        j2000[k,0] = datetime2j2000(dt)
    return j2000

def datetimearr2sod(datetimearr):
    """
    Converts a n x 1 array of python datetimes
    to an n x 1 second of day array
    """
    #Input sanitize
    if isinstance(datetimearr,numpy.ndarray):
        datetimearr = datetimearr.flatten().tolist() #Just easier to deal with lists

    sod = numpy.empty((len(datetimearr),1))
    for k,dt in enumerate(datetimearr):
        sod[k,0] = datetime2sod(dt)
    return sod

def datearr2sod(datearr):
    """
    Converts a n x 6 array with columns year,month,day,hour,minute,second
    to an n x 1 second of day array
    This and the next method show two ways of dealing with casting
    the input to int. Which is faster?
    """
    #Input sanitize
    if not isinstance(datearr,numpy.ndarray):
        raise RuntimeError('Input must be numpy.ndarray n x 6')

    sod = numpy.empty((len(datearr[:,0]),1))
    for k in range(len(datearr[:,0])):
        sod[k,0] = ymdhms2sod(int(datearr[k,0]),int(datearr[k,1]),int(datearr[k,2]),int(datearr[k,3]),int(datearr[k,4]),int(datearr[k,5]))
    return sod

def sodarr2datetime(sodarr,year,month,day):
    """
    Converts a nx1 array of seconds of the day
    to a datetime array. Year, month, day are
    scalar integers. Returns a numpy array of
    datetime objects.
    """

    #Input sanitize
    if isinstance(sodarr,numpy.ndarray):
        sodarr = sodarr.flatten().tolist() #Easier to deal with list, no indexing ambiguity

    dt = numpy.empty((len(sodarr),1),dtype='object')
    for k in range(len(sodarr)):
        dt[k,0] = sod2datetime(sodarr[k],year,month,day)
    return dt

def datearr2doy(datearr):
    """
    Converts a n x 6 array with columns year,month,day,hour,minute,second
    to an n x 1 decimal day of year array
    """

    #Input sanitize
    if not isinstance(datearr,numpy.ndarray):
        raise RuntimeError('Input must be numpy.ndarray n x 6')

    datearr = datearr.astype('int')
    doy = numpy.array((len(datearr[:,0]),1))
    for k in range(len(datearr[:,0])):
        doy[k,0] = ymdhms2doy(datearr[k,0],datearr[k,1],datearr[k,2],datearr[k,3],datearr[k,4],datearr[k,5])
    return doy

def datearr2datetime(datearr,asArray=True):
    """
    Converts a n x 6 array with columns year,month,day,hour,minute,second
    to an n x 1 list, or if asArray==True, as a numpy array
    This is sort of slow, but it gets the job done
    """
    datetimes = []
    datearr = datearr.astype('int')
    for k in range(len(datearr[:,0])):
        datetimes.append(datetime.datetime(datearr[k,0],datearr[k,1],datearr[k,2],datearr[k,3],datearr[k,4],datearr[k,5]))
    if asArray:
        datetimes = numpy.array(datetimes,ndmin=2)

    return datetimes

def datenumarr2datetime(datenumarr,asArray=True):
    """
    Converts a n x 1 array of matlab datenumbers to
    to an n x 1 list, or if asArray==True, as a numpy array
    This is sort of slow, but it gets the job done
    """

    datetimes = [datenum2datetime(datenumarr[k]) for k in range(len(datenumarr))]
    if asArray:
        datetimes = numpy.array(datetimes)
    return datetimes

def matchTimes(primary_dt,dt,tol_s=1,tol_us=4e5,fail_on_duplicates=True,allow_duplicates=False,warn_no_match=False):
    """
    Finds a matching timestamp in primary_dt
    within tolerance tol_us (given in microseconds)
    for every value in dt.
    Inputs:
    -------
        dt - numpy.array(dtype=object)
        primary_dt - numpy.array(dtype=object)

    Returns:
    --------
        inds - numpy.ndarray
        Array of len(dt) of indices into
        primary_dt or NaN if no match was found
        for a partiuclar dt value.
    """

    import bisect

    #Formatting String
    fmstr = '%Y-%m-%d %H:%M:%S.%f'

    # Input sanitize
    if not isinstance(dt,numpy.ndarray):
        dt = numpy.array(dt)
    if not isinstance(primary_dt,numpy.ndarray):
        primary_dt = numpy.array(primary_dt)

    # Throw an error if there are NaN values
    #if any(numpy.isnan(dt)) or any(numpy.isnan(primary_dt)):
    #   raise RuntimeError('NaN datetimes found! Sanitize your inputs before calling matchTimes!')
    #Nope...apparently isnan doesn't work on datetime arrays

    #Reshape things to be column vector if nessecary
    primary_dt=lmk_utils.asColumn(primary_dt)
    dt = lmk_utils.asColumn(dt)

    #Preallocate the results array and fill it with nan
    inds = lmk_utils.nan( (len(dt),1) )

    n_unmatched = 0
    for i in range(len(dt)):
        try:
            #Find the place where dt[i,0] should be inserted into primary_dt to maintain
            #sorted order. bisect and bisect_left should always return the same index
            #UNLESS the exact value of dt[i,0] is present
            ind = bisect.bisect(primary_dt[:,0].tolist(),dt[i,0])
            ind_l = bisect.bisect_left(primary_dt[:,0].tolist(),dt[i,0])

            if ind == None or ind_l == None:
                raise RuntimeError('Bisect returned None while trying to match %s!' % (dt[i,0].strftime(fmstr)))

            elif ind != ind_l: #This means that the exact value of dt[i] is present at ind_l
                inds[i] = ind_l

            elif ind > 0 and ind < len(primary_dt): #this means that primary_dt[ind] > dt[i] and primary_dt[ind-1] < dt[i]
                #Compute deltas and compare to tolerence
                delta_before = dt[i,0]-primary_dt[ind-1,0]
                delta_after = primary_dt[ind,0]-dt[i,0]

                #Compute before and after time difference in microseconds
                delta_before_us = abs(delta_before.total_seconds())*1.0e6
                delta_after_us = abs(delta_after.total_seconds())*1.0e6

                #Determine if the time difference to the index before/after is within tolerence
                before_intol = True if delta_before_us < tol_us else False
                after_intol = True if delta_after_us < tol_us else False

                #Determine if either of the values we'd like to assign is already matched
                before_duplicate = True if ind-1 in inds else False
                after_duplicate = True if ind in inds else False

                #Handle case of ambiguous match
                #Default behaviour is to choose best match if Ambigious
                if after_intol and before_intol:

                    #Determine which value provides the smallest difference in timestamp
                    best_match = 'before' if delta_before_us<=delta_after_us else 'after'

                    #Believe it or not, I've had the case happen where delta_before exactly equals delta_after. In
                    #this case, we assume that we want to use the earlier timestamp because of 'the arrow of time'
                    #i.e. using information from the past is less uncertain that from the future, all other things
                    #being equal

                    log.warn('Ambiguous match between two timestamps.\n%s within tolerance %f s of both\n %s (%f s) and\n %s (%f s)' % (dt[i,0].strftime(fmstr),tol_us/1.0e6,\
                            primary_dt[ind,0].strftime(fmstr),abs(delta_after.total_seconds()),
                            primary_dt[ind-1,0].strftime(fmstr),abs(delta_before.total_seconds())) )
                    if best_match=='before' and (not before_duplicate or allow_duplicates):
                        inds[i] = ind-1
                    elif best_match=='after' and (not after_duplicate or allow_duplicates):
                        inds[i] = ind
                else:
                    #No need to use best_match here, we only have one value within tolerance
                    if before_intol and (not before_duplicate or allow_duplicates):
                        inds[i] = ind-1
                    elif after_intol and (not after_duplicate or allow_duplicates):
                        inds[i] = ind
                    else:
                        if warn_no_match:
                            log.warn('No match found for timestamp %s' % (dt[i,0].strftime(fmstr)))
                        n_unmatched+=1
                        #log.warn('No match found for timestamp %s' % (dt[i,0].strftime(fmstr)))


            elif ind == 0 or ind == len(primary_dt):
                #Edge cases, use bisect left result if within tolerence
                if ind==len(primary_dt):
                    ind-=1
                delta = dt[i,0] - primary_dt[ind,0]
                if abs(delta.total_seconds()*1.0e6) < tol_us and (ind not in inds or allow_duplicates):
                    inds[i] = ind
                else:
                    n_unmatched+=1
                    #log.warn('No match found for timestamp %s [NOTE: Edge case, closest match was beginning or end of primary_dt]' % \
                    #   (dt[i,0].strftime(fmstr)))

            else:
                log.error('Unhandled bisect case...this should not happen')
                raise RuntimeError('Something wierd is going on, unhandled bisect case')

        except:
            traceback.print_exc()
            print("Error in matchTimes:\nIndex in dt: %d\nBisect Right Index in primary_dt: %d\nBisect Left Index in primary_dt: %d" % (i,ind,ind_l))
            print("dt[i] %s\nprimary_dt[ind] %s\nprimary_dt[ind_l] %s" % (str(dt[i]),str(primary_dt[ind]),str(primary_dt[ind_l])))
            pdb.set_trace()

    #print "Completed timestamp alignment, %d points remained unmatched out of %d\n" %(unmatched_counter,len(dt))

    inds = inds.astype('int')

    #Identify the first occurance of each unique value in inds in unqinds
    unqvals,unqinds = numpy.unique(inds,return_index=True)
    #Identifiy all indices in inds which are nan
    naninds = numpy.nonzero(numpy.isnan(inds))

    nunq = len(unqvals)
    log.info("%d/%d values to match were matched uniquely." % (nunq,len(inds)))
    log.info("%d/%d values are nan (unmatched)." % (len(naninds),len(inds)))

    return inds

def fastMatchTimes(primary_dt,dt,tol_us=4e5,fail_on_duplicates=True,allow_duplicates=False):
    """
    Finds a matching timestamp in primary_dt
    within tolerance tol_us (given in microseconds)
    for every value in dt. Faster implementation than above, but also
    less clear.
    Inputs:
    -------
        dt - numpy.array(dtype=object)
        primary_dt - numpy.array(dtype=object)

    Returns:
    --------
        inds - numpy.ndarray
        Array of len(dt) of indices into
        primary_dt or NaN if no match was found
        for a partiuclar dt value.
    """

    import bisect

    #Formatting String
    fmstr = '%Y-%m-%d %H:%M:%S.%f'

    reftime = datetime.datetime.utcfromtimestamp(0)

    # Input sanitize
    if isinstance(dt,numpy.ndarray):
        dt = dt.flatten().tolist()
    if isinstance(primary_dt,numpy.ndarray):
        primary_dt = primary_dt.flatten().tolist()

    #Convert to unix epoch (seconds since hour 0 Jan 1 1970)
    log.debug("Computing UNIX epoch milisecond integer values for all datetimes to improve matching speed")
    dt = numpy.array([(d-reftime).total_seconds() for d in dt])
    primary_dt = numpy.array([(d-reftime).total_seconds() for d in primary_dt])

    # Throw an error if there are NaN values
    #if any(numpy.isnan(dt)) or any(numpy.isnan(primary_dt)):
    #   raise RuntimeError('NaN datetimes found! Sanitize your inputs before calling matchTimes!')
    #Nope...apparently isnan doesn't work on datetime arrays

    #Preallocate the results array and fill it with nan
    inds = numpy.zeros((len(dt),1))
    inds.fill(numpy.nan)
    #inds = inds.flatten().tolist()

    starttime = datetime.datetime.now()

    ep2dt = lambda ep: reftime + datetime.timedelta(seconds=ep)

    log.debug("Beginning matching")
    n_unmatched = 0
    for i in range(len(dt)):
        if numpy.mod(i,5000)==0 and i>0:
            rate = float(i)/(datetime.datetime.now()-starttime).total_seconds()
            timetoend = (len(dt)-i)/rate/60.
            log.info("%d/%d matched, %d unmatched: ~%.3f iterations/sec. About %.1f minutes remain." % (i,len(dt),n_unmatched,rate,timetoend))
            n_unmatched = 0
        try:
            #Find the place where dt[i,0] should be inserted into primary_dt to maintain
            #sorted order. bisect and bisect_left should always return the same index
            #UNLESS the exact value of dt[i,0] is present
            ind = numpy.searchsorted(primary_dt,dt[i],side='right')

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
            pdb.set_trace()

    #print "Completed timestamp alignment, %d points remained unmatched out of %d\n" %(unmatched_counter,len(dt))

    #Convert list back to array
    #inds = numpy.array(inds)
    inds = inds.astype(int)
    #inds = lmk_utils.asColumn(inds)

    #Identify the first occurance of each unique value in inds in unqinds
    unqvals,unqinds = numpy.unique(inds,return_index=True)
    #Identifiy all indices in inds which are nan
    naninds = numpy.nonzero(numpy.isnan(inds))

    nunq = len(unqvals)
    log.info("%d/%d values to match were matched uniquely." % (nunq,len(inds)))
    log.info("%d/%d values are nan (unmatched)." % (len(naninds),len(inds)))

    return inds
