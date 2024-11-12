__all__ = [
    'plumerise_briggs', 'open_date', 'gd2hemco', 'pt2hemco', 'pt2gd', 'merge',
    'to_ioapi'
]

import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)

_defez = np.array([
    -6.00, 123.0, 254.0, 387.0, 521.0, 657.0, 795.0, 934.0, 1075., 1218.,
    1363., 1510., 1659., 1860., 2118., 2382., 2654., 2932., 3219., 3665.,
    4132., 4623., 5142., 5692., 6277., 6905., 7582., 8320., 9409., 10504,
    11578, 12633, 13674, 14706, 15731, 16753, 17773, 18807, 19855, 20920,
    22004, 23108, 24240, 25402, 26596, 27824, 29085, 30382, 31716, 33101,
    34539, 36030, 37574, 39173, 40825, 42529, 44286, 46092, 47946, 49844,
    51788, 53773, 55794, 57846, 59924, 62021, 64129, 66245, 68392, 70657,
    73180, 76357, 80581
], dtype='f')
_deflevs = np.array([
    0.99250002413, 0.97749990013, 0.962499776, 0.947499955, 0.93250006,
    0.91749991, 0.90249991, 0.88749996, 0.87249996, 0.85750006, 0.842500125,
    0.82750016, 0.8100002, 0.78750002, 0.762499965, 0.737500105, 0.7125001,
    0.6875001, 0.65625015, 0.6187502, 0.58125015, 0.5437501, 0.5062501,
    0.4687501, 0.4312501, 0.3937501, 0.3562501, 0.31279158, 0.26647905,
    0.2265135325, 0.192541016587707, 0.163661504087706, 0.139115, 0.11825,
    0.10051436, 0.085439015, 0.07255786, 0.06149566, 0.05201591, 0.04390966,
    0.03699271, 0.03108891, 0.02604911, 0.021761005, 0.01812435, 0.01505025,
    0.01246015, 0.010284921, 0.008456392, 0.0069183215, 0.005631801,
    0.004561686, 0.003676501, 0.002948321, 0.0023525905, 0.00186788,
    0.00147565, 0.001159975, 0.00090728705, 0.0007059566, 0.0005462926,
    0.0004204236, 0.0003217836, 0.00024493755, 0.000185422, 0.000139599,
    0.00010452401, 7.7672515e-05, 5.679251e-05, 4.0142505e-05, 2.635e-05,
    1.5e-05
], dtype='d')


_levattrs = dict(
    long_name="hybrid level at midpoints ((A/P0)+B)",
    units="level",
    axis="Z",
    positive="up",
    standard_name="atmosphere_hybrid_sigma_pressure_coordinate",
    formula_terms="a: hyam b: hybm p0: P0 ps: PS",
)
_latattrs = dict(
    long_name="Latitude",
    units="degrees_north",
    axis="Y",
    bounds="lat_bnds",
)
_lonattrs = dict(
    long_name="Longitude",
    units="degrees_east",
    axis="X",
    bounds="lon_bnds",
)
_reftime = '1970-01-01 00:00:00'
_timeattrs = dict(
    long_name="Time",
    units=f"hours since {_reftime}",
    calendar="gregorian",
    axis="T",
)


def getmw(key, gc='cb6r5_ae7_aq', nr='cb6r5hap_ae7_aq'):
    import requests
    import os
    import io
    import pandas as pd
    import re

    mwpath = f'cmaq_{gc}_molwt.csv'
    fillin = {
        'CH4': 16.0,         # from ECH4
        'ETHYLBENZ': 106.2,  # from XYLMN
        'BENZ': 78.1,        # from BENZENE
        'NH3': 17.0,         # from NR_{mech}.nml
    }
    if not os.path.exists(mwpath):
        mwdfs = []
        for prfx, mech in [('GC', gc), ('NR', nr)]:
            url = (
                'https://raw.githubusercontent.com/USEPA/CMAQ/refs/heads/main/'
                f'CCTM/src/MECHS/{mech}/{prfx}_{mech}.nml'
            )
            r = requests.get(url)
            txtlines = r.text.split('\n')
            datlines = [
                _l for _l in txtlines
                if _l.startswith("'") or _l.startswith('!')
            ]
            datlines[0] = datlines[0].replace('!', '') + ','
            dat = '\n'.join(datlines).replace("'", "")
            dat = re.sub('[ ]+,', ',', dat)
            rmwdf = pd.read_csv(io.StringIO(dat), index_col=False)
            rmwdf.columns = [k for k in rmwdf.columns]
            mwdfs.append(rmwdf.set_index('SPECIES'))
        mwdf = pd.concat(mwdfs)
        for newk, mw in fillin.items():
            if newk not in mwdf.index:
                mwdf.loc[newk, 'MOLWT'] = mw
        mwdf[['MOLWT']].to_csv(mwpath, index=True)
    mwdf = pd.read_csv(mwpath, index_col=0)
    try:
        mw = mwdf.loc[key, 'MOLWT'] / 1e3
    except KeyError:
        raise KeyError(f'{key} not found in {mwpath}')
    return mw


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


def pt2hemco(
    path, pf, elat, elon, ez=None, nk=11, temp_a=288.15, pres_a=101325, u=2.5
):
    """
    Arguments
    ---------
    path : str
        Path to save as HEMCO file
    pf : xarray.Dataset
        Point source from se_file
    elat : array
        Edge latitudes for regular grid
    elon : array
        Edge longitudes for regular grid
    ez : array
        Edge altitudes in meeters for vertical structure
    nk : int
        Number of vertical levels to use if ek not specified
    temp_a : float
        Temperature in K for temperature to use for plume rise.
    pres_a : float
        Pressure in Pa for temperature to use for plume rise.

    Arguments
    ---------
    outf : hemcofile
        Object with .nc property
    """
    import numpy as np
    import pandas as pd
    assert pf['lat'].min() > elat.min()
    assert pf['lat'].max() < elat.max()
    assert pf['lon'].min() > elon.min()
    assert pf['lon'].max() < elon.max()

    # allocating emissions to level based on stack height (STKHT) and Briggs
    # plume rise using approximate level heights from
    # http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
    # and assuming constant T=288.15 and P=101325
    if ez is None:
        ez = _defez[:nk + 1]
    nt = pf.sizes['time']
    nr = elat.size - 1
    nc = elon.size - 1
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    ilat = np.arange(nr)
    ilon = np.arange(nc)
    dx = pf.attrs['XCELL'] / 2
    # replace continues lat/lon with midpoint indexer
    ris = pd.cut(pf['lat'], bins=elat, labels=ilat).astype('i')
    cis = pd.cut(pf['lon'], bins=elon, labels=ilon).astype('i')
    if ez is not None and 'STKHT' in pf:
        # cz = (ez[1:] + ez[:-1]) / 2
        iz = np.arange(nk)
        dz = plumerise_briggs(
            pf['STKDM'], pf['STKVE'], pf['STKTK'],
            temp_a=temp_a, pres_a=pres_a, u=u, x=dx
        )
        z = pf['STKHT'] + dz
        z = np.minimum(np.maximum(z, ez[0]), ez[-1])
        kis = pd.cut(z, bins=ez, labels=iz).astype('i')
    else:
        nk = 1
        kis = np.zeros(pf.sizes['stack'], dtype='i')

    clev = _deflevs[:nk]
    tis = (
        (pf.time - pf.time.min()).dt.total_seconds() / 3600
    ).round(0).astype('i').data
    pf['ti'] = ('time',), tis
    pf['ki'] = ('stack',), kis
    pf['ri'] = ('stack',), ris
    pf['ci'] = ('stack',), cis
    nt = 25
    tmp = np.zeros((nt, nk, nr, nc), dtype='f')
    datakeys = [
        k for k, v in pf.data_vars.items()
        if (
            k not in ('TFLAG', 'lon', 'lat', 'ti', 'ki', 'ri', 'ci')
            and len(v.dims) > 1
        )
    ]
    outf = hemcofile(
        path, pf.time, clat, clon, lev=clev, varkeys=datakeys, attrs=pf.attrs
    )
    area = hemco_area(elat, elon)
    outf.addvar('AREA', area, units='m2', dims=('lat', 'lon'))
    for dk in datakeys:
        if len(pf[dk].dims) == 1:
            print(f'skip {dk}')
            continue
        dtot = pf[dk].sum()
        if dtot == 0:
            print(f'excluded {dk}')
            continue
        tmp[:] *= 0
        print(dk)
        df = pf[['ti', 'ki', 'ri', 'ci', dk]].to_dataframe()
        df = df.loc[df[dk] != 0]
        vals = df.groupby(['ti', 'ki', 'ri', 'ci'], as_index=False).sum()
        tmp[vals.ti, vals.ki, vals.ri, vals.ci] = vals[dk].values
        attrs = {k: v for k, v in pf[dk].attrs.items()}
        unit = attrs['units'].strip()
        tmp, unit = unitconvert(dk, tmp, unit, area=area)
        attrs['units'] = unit
        outf.addvar(dk, tmp, **attrs)

    return outf


def gd2matrix(gf, elat, elon):
    from shapely.geometry import box
    import geopandas as gpd
    import warnings
    clat = (elat[:-1] + elat[1:]) / 2
    clon = (elon[:-1] + elon[1:]) / 2
    hdx = np.diff(elon).mean() / 2
    hdy = np.diff(elat).mean() / 2
    qgeodf = gf.csp.geodf
    qgeodf['original_area'] = qgeodf.geometry.area
    if hdx == hdy:
        latj = np.arange(clat.size)
        loni = np.arange(clon.size)
        LONI, LATJ = np.meshgrid(loni, latj)
        LON, LAT = np.meshgrid(clon, clat)
        LAT = LAT.ravel()
        LON = LON.ravel()
        LATJ = LATJ.ravel()
        LONI = LONI.ravel()
        hgeodf = gpd.GeoDataFrame(
            dict(lat=LAT, lon=LON, lati=LATJ, loni=LONI),
            geometry=gpd.points_from_xy(LON, LAT), crs=4326
        )
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            hgeodf['geometry'] = hgeodf['geometry'].buffer(
                hdx, cap_style='square'
            )
    else:
        geoms = []
        lats = []
        latis = []
        lons = []
        lonis = []
        for lati, lat in enumerate(clat):
            for loni, lon in enumerate(clon):
                lats.append(lat)
                latis.append(lati)
                lons.append(lons)
                lonis.append(loni)
                geoms.append(box(lon - hdx, lat - hdy, lon + hdx, lat + hdy))
        hgeodf = gpd.GeoDataFrame(
            dict(lat=lats, lons=lons, lati=latis, loni=lonis),
            geometry=geoms, crs=4326
        )
    ol = gpd.overlay(qgeodf.reset_index(), hgeodf.to_crs(qgeodf.crs))
    ol['intx_area'] = ol.geometry.area
    ol['fraction'] = ol['intx_area'] / ol['original_area']
    ol['ri'] = ol['ROW'].astype('i')
    ol['ci'] = ol['COL'].astype('i')
    return ol.set_index(['ROW', 'COL', 'lati', 'loni'])


def gd2hemco_fast(path, gf, elat, elon):
    """
    Bilinear interpolation of fluxes (w/ MSFX2 factor)

    Arguments
    ---------
    path : str
        Path to save as HEMCO file
    pf : xarray.Dataset
        Point source from se_file
    elat : array
        Edge latitudes for regular grid
    elon : array
        Edge longitudes for regular grid

    Arguments
    ---------
    outf : hemcofile
        Object with .nc property
    """
    import pyproj
    proj = pyproj.Proj(gf.crs)
    qarea = (
        gf.XCELL * gf.YCELL
        / proj.get_factors(gf['lon'], gf['lat']).areal_scale
    )
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    clev = _deflevs[:1]
    datakeys = [
        k for k in gf
        if k not in ('TFLAG', 'ti', 'ki', 'ri', 'ci', 'lon', 'lat')
    ]
    outf = hemcofile(
        path, gf.time, clat, clon, lev=clev, varkeys=datakeys, attrs=gf.attrs
    )
    area = hemco_area(elat, elon)
    outf.addvar('AREA', area, units='m2', dims=('lat', 'lon'))
    LON, LAT = np.meshgrid(clon, clat)
    X, Y = proj(LON, LAT)
    X = xr.DataArray(X, dims=('ROW', 'COL'))
    Y = xr.DataArray(Y, dims=('ROW', 'COL'))
    for dk in datakeys:
        dtot = gf[dk].sum()
        if dtot == 0:
            print(f'excluded {dk}')
            continue
        print(dk)
        tmp = (gf[dk] / qarea).interp(ROW=Y, COL=X)
        attrs = {k: v for k, v in gf[dk].attrs.items()}
        unit = attrs['units'].strip() + '/m**2'
        tmp, unit = unitconvert(dk, tmp, unit=unit)
        attrs['units'] = unit
        outf.addvar(dk, tmp.data, **attrs)

    return outf


def gd2hemco(path, gf, elat, elon, matrix=None):
    """
    Uses a fractional aera overlap interoplation.

    Arguments
    ---------
    path : str
        Path to save as HEMCO file
    pf : xarray.Dataset
        Point source from se_file
    elat : array
        Edge latitudes for regular grid
    elon : array
        Edge longitudes for regular grid
    matrix : pandas.DataFrame
        fraction from row/col centroids to lat/lon centroids

    Arguments
    ---------
    outf : hemcofile
        Object with .nc property
    """
    import numpy as np
    import pandas as pd
    if not gf['lat'].min() > elat.min():
        raise AssertionError('input grid lat lower than output grid')
    if not gf['lat'].max() < elat.max():
        raise AssertionError('input grid lat greater than output grid')
    if not gf['lon'].min() > elon.min():
        raise AssertionError('input grid lon lower than output grid')
    if not gf['lon'].max() < elon.max():
        print(gf['lon'].max())
        raise AssertionError('input grid lon greater than output grid')
    if matrix is None:
        matrix = gd2matrix(gf, elat, elon)

    nt = gf.sizes['time']
    nr = elat.size - 1
    nc = elon.size - 1
    nk = 1
    clat = (elat[1:] + elat[:-1]) / 2
    clon = (elon[1:] + elon[:-1]) / 2
    clev = _deflevs[:nk]
    tis = (
        (gf.time - gf.time.min()).dt.total_seconds() / 3600
    ).round(0).astype('i').data
    gf['ti'] = ('time',), tis
    kis = np.zeros((nk,), dtype='i')
    gf['ki'] = ('LAY',), kis
    nt = 25
    tmp = np.zeros((nt, nk, nr, nc), dtype='f')
    datakeys = [
        k for k in gf
        if k not in ('TFLAG', 'ti', 'ki', 'ri', 'ci', 'lon', 'lat')
    ]
    outf = hemcofile(
        path, gf.time, clat, clon, lev=clev, varkeys=datakeys, attrs=gf.attrs
    )
    area = hemco_area(elat, elon)
    outf.addvar('AREA', area, units='m2', dims=('lat', 'lon'))
    if matrix is not None:
        matrix = matrix[['fraction']].reset_index()
    else:
        matrix = gf[['ROW', 'COL']].to_dataframe().reset_index()
        ilat = np.arange(nr + 1)
        ilon = np.arange(nc + 1)
        # replace continues lat/lon with midpoint indexer
        ris = np.interp(
            gf['lat'], elat, ilat, left=-999, right=-888
        ).astype('i')
        cis = np.interp(
            gf['lon'], elon, ilon, left=-999, right=-888
        ).astype('i')
        matrix['lati'] = ris.data.ravel()
        matrix['loni'] = cis.data.ravel()
        matrix['fraction'] = 1

    merged = None
    for dk in datakeys:
        if len(gf[dk].dims) == 1:
            print(f'skip {dk}')
            continue
        tmp[:] *= 0
        print(dk)
        df = gf[['ti', 'ki', dk]].to_dataframe()
        gvals = df.groupby(['ti', 'ki', 'ROW', 'COL'], as_index=True).sum()
        if merged is None:
            merged = pd.merge(
                gvals.reset_index(), matrix,
                left_on=['ROW', 'COL'], right_on=['ROW', 'COL'], how='left'
            ).rename(columns={dk: 'val'}).set_index(
                ['ti', 'ki', 'ROW', 'COL', 'lati', 'loni']
            )
            midx = merged.index.droplevel(['lati', 'loni'])
            frac = merged['fraction'].values
        merged['val'] = frac * gvals.loc[midx, dk].values
        idxs = ['ti', 'ki', 'lati', 'loni']
        gval = merged['val'].groupby(idxs).sum()
        ti = gval.index.get_level_values('ti')
        ki = gval.index.get_level_values('ki')
        lati = gval.index.get_level_values('lati')
        loni = gval.index.get_level_values('loni')
        tmp[ti, ki, lati, loni] += gval
        attrs = {k: v for k, v in gf[dk].attrs.items()}
        tmp, unit = unitconvert(dk, tmp, attrs['units'], area=area)
        attrs['units'] = unit
        outf.addvar(dk, tmp, **attrs)

    return outf


def unitconvert(key, val, unit, area=None, inplace=True):
    outunit = []
    factor = np.ones_like(val)
    assert '/s' in unit
    inunit = unit.replace('/s', '')
    if unit in ('g/s', 'g/s/m**2'):
        factor /= 1000.
        outunit.append('kg')
    elif unit in ('moles/s', 'moles/s/m**2'):
        try:
            mw = getmw(key)
            factor *= mw
            outunit.append('kg')
        except KeyError as e:
            print(f'**WARNING: {key} in {unit} not converted to kg: {e}')
            outunit.append(inunit)
    else:
        print(f'**WARNING: {key} [{unit}] not converted to kg: {unit} unknown')
        outunit.append(inunit)
    if 'm**2' not in unit and area is not None:
        factor /= area
        outunit.append('/m**2')
    elif '/m**2' in unit:
        outunit.append('/m**2')
    outunit.append('/s')
    outunit = ''.join(outunit)
    if inplace:
        outval = val
        outval *= factor
    else:
        outval = val * factor
    outunit = outunit.replace('/s/m', '/m2/s')
    return outval, outunit


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
    """
    Convert point file (from se_file) to a CMAQ gridded emission file.
    """
    import numpy as np
    import xarray as xr
    import pandas as pd

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


def hemco_area(elat, elon, R=6371007.2):
    """
    Arguments
    ---------
    elat : array
        Edges that define a regular grid in degrees from 0N
    elon : array
        Edges that define a regular grid in degrees from 0E
    R : float
        Radius of the earth https://github.com/geoschem/geos-chem/issues/2453

    Returns
    -------
    A : array
        Area in square meters with dimensions (nlat, nlon)
    """
    import numpy as np
    dx = np.diff(elon).mean()
    # 1-d area
    a = dx * np.pi / 180 * R**2 * np.diff(np.sin(np.radians(elat)))
    # 2-d area
    A = a.astype('f')[:, None].repeat(elon.size - 1, 1)
    return A


class hemcofile:
    def __init__(
        self, path, time, lat, lon, lev=None, varkeys=None, attrs=None
    ):
        import pandas as pd
        import netCDF4
        if varkeys is None:
            varkeys = []
        nc = self.nc = netCDF4.Dataset(
            path, mode='ws', format='NETCDF4_CLASSIC'
        )
        if attrs is not None:
            for pk, pv in attrs.items():
                nc.setncattr(pk, pv)
        nc.createDimension('time', None)
        if lev is not None:
            nc.createDimension('lev', lev.size)
        nc.createDimension('lat', lat.size)
        nc.createDimension('lon', lon.size)
        tv = nc.createVariable('time', 'd', ('time',))
        for k, v in _timeattrs.items():
            tv.setncattr(k, v)
        timec = (
            (
                pd.to_datetime(time)
                - pd.to_datetime(_reftime)
            ).total_seconds() / 3600.
        ).astype('d')
        tv[:timec.size] = timec
        if lev is not None:
            levv = nc.createVariable('lev', 'd', ('lev',))
            for k, v in _levattrs.items():
                levv.setncattr(k, v)
            levv[:] = lev
        latv = nc.createVariable('lat', 'd', ('lat',))
        for k, v in _latattrs.items():
            latv.setncattr(k, v)
        latv[:] = lat
        lonv = nc.createVariable('lon', 'd', ('lon',))
        for k, v in _lonattrs.items():
            lonv.setncattr(k, v)
        lonv[:] = lon
        for vk in varkeys:
            self.defvar(vk)

    def defvar(self, vk, dims=None, **attrs):
        """
        Define a variable using HEMCO expectations
        """
        ncf = self.nc
        nr = len(ncf.dimensions['lat'])
        nc = len(ncf.dimensions['lon'])
        chunkdefsizes = dict(time=1, lev=1, lat=nr, lon=nc)
        if dims is None:
            if 'lev' not in ncf.dimensions:
                dims = ('time', 'lat', 'lon')
            else:
                dims = ('time', 'lev', 'lat', 'lon')
        chunks = [chunkdefsizes[dk] for dk in dims]
        vv = ncf.createVariable(
            vk, 'f', dims, chunksizes=chunks, zlib=True, complevel=1
        )
        vv.setncattr('standard_name', vk)
        vv.setncattr('units', 'unknown')
        for pk, pv in attrs.items():
            vv.setncattr(pk, pv)

    def setattrs(self, **attrs):
        """
        Add attributes to file
        """
        for k, v in attrs.items():
            self.nc.setncattr(k, v)

    def addvar(self, key, vals, dims=None, **attrs):
        """
        Add a variable (defining if necessary using HEMCO expectations)
        """
        nc = self.nc
        if key not in nc.variables:
            self.defvar(key, dims=dims, **attrs)
        vv = nc.variables[key]
        for pk, pv in attrs.items():
            vv.setncattr(pk, pv)
        nt = vals.shape[0]
        vv[:nt] = vals

    def close(self):
        self.nc.close()

    def __del__(self):
        self.nc.close()
        del self.nc


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


def symlinks(tmpl, dates, datetype=None, verbose=0):
    """
    Arguments
    ---------
    tmpl : str
        strftime template for paths
    dates : pd.Series or str
        If Series, must have date index (destination date) and date values
        (source)
        If str, path to csv file with dates file as describe in Notes.
    datetype: str or None
        If dates is an instance of str, choose one of: aveday_Y, aveday_N,
        mwdss_Y, mwdss_N, week_Y, week_N, all. Where Y/N denotes if "holidays"
        have special treatment.
    verbose: int
        Level of verbosity

    Returns
    -------
    links : list
        List of links that were created

    Notes
    -----
    Date,aveday_N,aveday_Y,mwdss_N,mwdss_Y,week_N,week_Y,all
    20160401,20160405,20160405,20160405,20160405,20160408,20160408,20160401
    ...

    """
    import pandas as pd
    import os
    if not isinstance(dates, pd.Series):
        datepath = dates
        df = pd.read_csv(datepath)
        df.columns = [k.strip() for k in df.columns]
        inkey = datetype
        outkey = df.columns[0]
        df[inkey] = pd.to_datetime(df[inkey])
        df[outkey] = pd.to_datetime(df[outkey])
        dates = df[[inkey, outkey]].set_index(outkey)
    links = []
    for outdate, indate in dates.items():
        src = indate.strftime(tmpl)
        dst = outdate.strftime(tmpl)
        if not os.path.exists(dst):
            if not os.path.exists(src):
                if verbose > 0:
                    print(f'{src} is missing; cannot make {dst}')
                continue
            os.symlink(src, dst)
            links.append(dst)
    return links
