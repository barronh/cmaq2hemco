__all__ = [
    'plumerise_briggs', 'open_date', 'grid2grid', 'pt2gd', 'merge', 'to_ioapi'
]

import xarray as xr
xr.set_options(keep_attrs=True)


def plumerise_briggs(
    stkdm, stkvel, stktk, pres_a=101325., temp_a=288.15, u=2.5, x=6000.,
    theta_lapse=None, F=None
):
    """
    Briggs (1969, 1971, 1974) equations of Plume Rise as documented within
    Seinfeld and Pandis[1]. On pg 868, Table 18.4 gives 7 equations -- 3 for
    stable conditions and 4 for neutral/unstable conditions. Stable
    calculations are only done with theta_lapse is provided.

    Arguments
    ---------
    stkdm : float
        Diameter of stack opening (m)
    stkvel : float
        Velocity of stack gas at opening (m/s)
    stktk : float
        Temperature of gas at opening (K)
    pres_a : float
        Pressure (Pa) of ambient environment; default 101325.
    temp_a : float
        Temperature (K) of ambient environment; default 288.13K
    u : float
        Wind speed (m/s); default of 2.5 m/s is used as a low estimate.
    x : float
        Distance from stack at which plume-rise is calculated (m). Default 6000
        is used because it is half of the commonly used grid cell size 12km.
    theta_lapse : float
        Potential Temperature Gradient (dtheta / dz) in K/m. Values below were
        converted from Seinfeld and Pandis (2006) Table 18.5 by dividing lapse
        rates per 100m by 100 to get the per meter to lapse rate by class.
          * A Extremely unstable dtheta / dz: < -0.009
          * B Moderately unstable dtheta / dz: -0.009 to -0.007
          * C Slightly unstable dtheta / dz: -0.007 to -0.005
          * D Neutral dtheta / dz: -0.005 to 0.005
          * E slightly stable dtheta / dz: 0.005 to 0.025
          * F moderately stable dtheta / dz : > 0.025
        If theta_lapse is none or less than or equal to 0.005, then the
        calculations for the unstable/neutral atmosphere are used. If
        theta_lapse is greater than 0.005, the three stable equations will be
        solved and the minimum will be returned.
    F : float
        Buoyancy parameter (m4/s3). If F is provided, it overrides the default
        calculation from stkdm, stkvel, and stktk.
    Returns
    -------
    dz : float
        Plume rise height of centerline.

    Notes
    -----
    Approximates calculations used by CMAQ and SMOKE to calculate plume rise.

    References
    ----------
    [1] Seinfeld, J. H. and Pandis, S. N.: Atmospheric chemistry and physics:
    from air pollution to climate change, 2nd ed., J. Wiley, Hoboken, N.J, 2006


    """
    import numpy as np
    g = 9.80665  # scipy.constants.g
    # Buoyancy flux parameter in Table 18.4 footnote
    if F is None:
        # do not allow negative temperature in F parameter
        dT = np.maximum(0, stktk - temp_a)
        F = g * stkdm**2 / 4 * stkvel * dT / stktk
    # As implemented on pg 868 of Seinfeld and Pandis ACP (2006)
    # using Neutral and unstable environments.
    if theta_lapse is None or theta_lapse <= 0.005:
        dzbriggs_neutral_loF_loX = (1.6 * F**(1 / 3.)) * x**(2 / 3.) / u
        dzbriggs_neutral_loF_hiX = (21.4 * F**(3 / 4.)) * x**(0) / u
        dzbriggs_neutral_hiF_loX = dzbriggs_neutral_loF_loX
        dzbriggs_neutral_hiF_hiX = (38.7 * F**(3 / 5.)) * x**(0) / u
        dzbriggs_loF = np.where(
            x < (49*F**(5/8.)),
            dzbriggs_neutral_loF_loX,
            dzbriggs_neutral_loF_hiX
        )
        dzbriggs_hiF = np.where(
            x < (119*F**(2/5.)),
            dzbriggs_neutral_hiF_loX,
            dzbriggs_neutral_hiF_hiX
        )
        dzbriggs = np.where(F < 55, dzbriggs_loF, dzbriggs_hiF)
    else:
        S2 = theta_lapse * g / temp_a  # theta_lapse
        dzbriggs_stable_1 = (2.4 * (F / S2))**(1 / 3.) / u**(1 / 3.)
        dzbriggs_stable_2 = (5.0 * F**(1 / 4.) * S2**(-3 / 8.))
        dzbriggs_stable_3 = 1.6 * F**(1 / 3.) * x**(2 / 3.) / u
        # The minimum of three equations is selected (Table 18.4 footnote c)
        dzbriggs = np.minimum(dzbriggs_stable_1, dzbriggs_stable_2)
        dzbriggs = np.minimum(dzbriggs, dzbriggs_stable_3)

    return dzbriggs


def open_date(
    date, tmpl, bucket, cache=True
):
    """
    Open all files for specific date

    Arguments
    ---------
    date : str
        Date parsable by pandas.to_datetime
    tmpl : str
        strftime template for date file
        (e.g., MCIP/12US1/GRIDCRO2D.12US1.35L.%y%m%d)
    bucket : str
        Bucket to pull from (e.g., )
    cache : bool
        Store file to prevent re-downloading.

    Returns
    -------
    ds : xarray.Dataset
        File opened (either in memory or from disk)
    """
    import boto3
    import pandas as pd
    from botocore import UNSIGNED
    from botocore.client import Config
    import io
    import os
    import gzip
    import cmaqsatproc as csp
    global res
    client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    date = pd.to_datetime(date)
    path = date.strftime(tmpl)
    dest = os.path.join(bucket, path)
    if cache:
        if not os.path.exists(dest):
            res = client.list_objects(Bucket=bucket, Prefix=path)
            nres = len(res.get('Contents', []))
            if nres != 1:
                raise IOError(f'{path} does not exist')
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            client.download_file(bucket, path, dest)

        if dest.endswith('.gz'):
            bdy = io.BytesIO(gzip.open(dest).read())
            f = csp.open_ioapi(bdy, engine='scipy')
        else:
            f = csp.open_ioapi(dest)
    else:
        res = client.list_objects(Bucket=bucket, Prefix=path)
        nres = len(res.get('Contents', []))
        if nres != 1:
            raise IOError(f'{path} does not exist')
        obj = client.get_object(Bucket=bucket, Key=path)
        bdy = obj['Body'].read()
        if path.endswith('.gz'):
            bdy = gzip.decompress()
        bdy = io.BytesIO(bdy)
        f = csp.open_ioapi(bdy, engine='scipy')
    return f


def grid2grid(f, elat, elon, ez=None):
    import numpy as np
    import pandas as pd
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    df = f.to_dataframe()
    # replace continues lat/lon with midpoint indexer
    df['lat'] = pd.cut(df['lat'], bins=elat, labels=clat)
    df['lon'] = pd.cut(df['lon'], bins=elon, labels=clon)
    if ez is not None:
        cz = (ez[1:] + ez[:-1]) / 2
        dz = plumerise_briggs(df['STKDM'], df['STKVE'], df['STKTK'])
        z = df['STKHT'] + dz
        z = np.minimum(np.maximum(z, ez[0]), ez[-1])
        df['lev'] = pd.cut(z, bins=ez, labels=cz)
        gbkeys = ['time', 'lev', 'lat', 'lon']
    else:
        gbkeys = ['time', 'lat', 'lon']
    # keep bins with no values and sum all emissions in each t/y/x bin.
    outdf = df.groupby(gbkeys, observed=False).sum(numeric_only=True)
    outds = outdf.to_xarray()
    chunks = [dl for dk, dl in outds.sizes.items()]
    chunks[0] = 1
    outds.attrs.update(f.attrs)
    for k in list(outds.data_vars):
        outds[k].encoding.update(complevel=1, chunks=chunks, zlib=True)
        outds[k].attrs.update(f[k].attrs)
    return outds


def merge(fs, bf=None):
    import copy
    if bf is None:
        bf = sorted([
            (f.sizes.get('lev', 1), f)
            for f in fs
        ])[-1][1][['time', 'lev', 'lat', 'lon']]
    fattrs = copy.deepcopy(bf.attrs)
    vattrs = {k: copy.deepcopy(v.attrs) for k, v in bf.data_vars.items()}
    for f in fs:
        for k in f:
            if k not in bf:
                bf[k] = f[k]
                vattrs[k] = copy.deepcopy(f[k].attrs)
            else:
                bf[k] = bf[k] + f[k].data
            bf[k].attrs.update(vattrs[k])
    bf.attrs.update(fattrs)
    return bf


def pt2gd(
    pf, nr, nc, ez=None, vgtyp=-9999, vgtop=5000., vglvls=None, byvar=True
):
    import numpy as np
    import xarray as xr
    import pandas as pd
    from .utils import plumerise_briggs

    ns = pf.sizes['stack']
    nt = pf.sizes['time']
    if ez is None:
        nz = 1
        kis = np.zeros((ns,), dtype='i')
    else:
        nz = ez.size - 1
        dz = plumerise_briggs(pf['STKDM'], pf['STKVE'], pf['STKTK'])
        dz[np.isnan(dz)] = 0
        zs = pf['STKHT'] + dz
        zs = np.minimum(ez.max(), np.maximum(ez.min(), zs))
        cz = np.arange(nz, dtype='i')
        kis = pd.cut(zs, bins=ez, labels=cz)
        kis = kis.astype('i')

    ex = np.arange(nc) * pf.XCELL + pf.XORIG
    ey = np.arange(nr) * pf.YCELL + pf.YORIG
    cx = np.arange(ex.size - 1, dtype='i')
    cy = np.arange(ey.size - 1, dtype='i')
    ris = pd.cut(pf['YLOCA'], bins=ey, labels=cy)
    cis = pd.cut(pf['XLOCA'], bins=ex, labels=cx)
    outside = (ris != ris) | (cis != cis)
    outf = xr.Dataset()
    outf.attrs.update(pf.attrs)
    outf.attrs['NCOLS'] = nc
    outf.attrs['NROWS'] = nr
    outf.attrs['NLAYS'] = nz
    if vglvls is None:
        vglvls = np.arange(ez.size)
    outf.attrs['VGLVLS'] = vglvls
    outf.attrs['VGTYP'] = vgtyp
    outf.attrs['VGTOP'] = float(vgtop)
    datakeys = [k for k in pf if k not in ('TFLAG',)]
    tmp = np.zeros((nt, nz, nr, nc), dtype='f')
    imod = max(1, ns // 1000)

    if byvar:
        pf['ti'] = ('time',), pf.time.dt.hour.data
        pf['ki'] = ('stack',), kis
        pf['ri'] = ('stack',), ris
        pf['ci'] = ('stack',), cis
    for dk in datakeys:
        if len(pf[dk].dims) == 1:
            print(f'\nskip {dk}')
            continue
        tmp[:] *= 0
        if byvar:
            print(dk)
            df = pf[['ti', 'ki', 'ri', 'ci', dk]].to_dataframe()
            df = df.loc[df[dk] != 0]
            vals = df.groupby(['ti', 'ki', 'ri', 'ci'], as_index=False).sum()
            tmp[vals.ti, vals.ki, vals.ri, vals.ci] = vals[dk].values
        else:
            vals = pf[dk].load()
            for si in np.arange(ns):
                if outside[si]:
                    print(f'\nskip {si}')
                    continue
                if (si % imod) == 0:
                    print(f'\r{dk:16s}: {si / ns:7.1%}', end='', flush=True)
                ki = int(kis[si])
                ri = int(ris[si])
                ci = int(cis[si])
                tmp[:, ki, ri, ci] += vals[:, si]
            print(f'\r{dk:16s}: {1:7.1%}', flush=True)

        outf[dk] = ('TSTEP', 'LAY', 'ROW', 'COL'), tmp, pf[dk].attrs
        outf[dk].encoding.update(pf[dk].encoding)
        outf[dk].encoding['chunks'] = (1, 1, nr, nc)

    return outf


def se_file(sf, ef):
    """
    Arguments
    ---------
    sf : xarray.Dataset
        Expected to be read from a CMAQ stack file (TSTEP, LAY, ROW, COL)
    ef : xarray.Dataset
        Expected to be read from a CMAQ emln file (TSTEP, LAY, ROW, COL)

    Returns
    -------
    cf : xarray.Dataset
        Combined file with emissions variables with dimensions ('time',
        'stack') and stack properties with dimensions ('stack',)
    """
    ef = ef.rename(ROW='stack', TSTEP='time')
    sf = sf.isel(TSTEP=0, LAY=0, COL=0, drop=True).rename(ROW='stack')
    del sf['TFLAG']
    del ef['TFLAG']
    for k in sf.data_vars:
        if k not in ef:
            ef[k] = sf[k]
    ef['lat'] = sf['LATITUDE']
    ef['lon'] = sf['LONGITUDE']
    del ef['LATITUDE']
    del ef['LONGITUDE']
    return ef


def to_ioapi(ef, path, **wopts):
    import xarray as xr
    import numpy as np
    import pandas as pd
    wopts.setdefault('mode', 'ws')
    wopts.setdefault('format', 'NETCDF4_CLASSIC')
    wopts.setdefault('unlimited_dims', ('TSTEP',))
    if 'TFLAG' not in ef:
        if 'time' in ef:
            time = ef.time.dt
        else:
            nt = ef.sizes['TSTEP']
            assert ef.attrs['TSTEP'] in (0, 10000)
            dt = pd.to_timedelta('1h') * np.arange(nt)
            time = pd.to_datetime([ef.SDATE] * nt, format='%Y%j') + dt
        date = time.strftime('%Y%j').astype('i')
        time = time.strftime('%H%M%S').astype('i')
        tflag = xr.DataArray(
            np.array([date, time]).T,  # t, 2
            dims=('TSTEP', 'DATE-TIME'),
            attrs=dict(
                units='<YYYYJJJ,HHMMSS>', long_name='TFLAG'.ljust(16),
                var_desc='TFLAG'.ljust(80)
            )
        ).expand_dims(VAR=np.arange(len(ef.data_vars))).transpose(
            'TSTEP', 'VAR', 'DATE-TIME'
        )
        ef['TFLAG'] = tflag
    ef.to_netcdf(path, **wopts)
