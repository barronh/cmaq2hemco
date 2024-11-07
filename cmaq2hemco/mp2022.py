__all__ = ['open_gdemis', 'open_ptemis']

import xarray as xr
from .utils import open_date
xr.set_options(keep_attrs=True)


def open_gdemis(date, sector):
    import pyproj
    eroot = 'emis/2022v1/2022hc_cb6_22m/12US1/cmaq_cb6ae7/premerged_area'
    bucket = 'epa-2022-modeling-platform'
    try:
        epath = (
            f'{eroot}/{sector}/emis_mole_{sector}_%Y%m%d'
            + '_12US1_cmaq_cb6ae7_2022hc_cb6_22m.ncf'
        )
        ef = open_date(date, epath, bucket).isel(
            LAY=0, drop=True
        ).rename(TSTEP='time')
    except Exception:
        epath = (
            f'{eroot}/{sector}/emis_mole_{sector}_%Y%m%d'
            + '_12US1_cmaq_cb6ae7_2022hc_cb6_22m.ncf.gz'
        )
        ef = open_date(date, epath, bucket).isel(
            LAY=0, drop=True
        ).rename(TSTEP='time')
    proj = pyproj.Proj(ef.crs)
    Y, X = xr.broadcast(ef.ROW, ef.COL)
    LON, LAT = proj(X, Y, inverse=True)
    attrs = dict(units='degrees_east', long_name='longitude')
    ef['lon'] = ('ROW', 'COL'), LON, attrs
    attrs = dict(units='degrees_north', long_name='latitude')
    ef['lat'] = ('ROW', 'COL'), LAT, attrs
    del ef['TFLAG']

    return ef


def open_ptemis(date, sector):
    from .utils import se_file
    eroot = 'emis/2022v1/2022hc_cb6_22m/12US1/cmaq_cb6ae7'
    bucket = 'epa-2022-modeling-platform'
    try:
        spath = (
            f'{eroot}/{sector}/stack_groups_{sector}_%Y%m%d'
            + '_12US1_2022hc_cb6_22m.ncf'
        )
        sf = open_date(date, spath, bucket)
    except Exception:
        spath = (
            f'{eroot}/{sector}/stack_groups_{sector}_12US1_2022hc_cb6_22m.ncf'
        )
        sf = open_date(date, spath, bucket)

    epath = (
        f'{eroot}/{sector}/inln_mole_{sector}_%Y%m%d'
        + '_12US1_cmaq_cb6ae7_2022hc_cb6_22m.ncf'
    )
    ef = open_date(date, epath, bucket).isel(
        LAY=0, COL=0, drop=True
    )
    ef = se_file(sf, ef)
    return ef


def process_dates(dates, gkeys, pkeys, elat, elon, ez, debug=False):
    import os
    from cmaq2hemco.utils import pt2gd, to_ioapi
    from cmaq2hemco.mp2022 import open_gdemis, open_ptemis
    # l1z = ez[:2].mean()
    nr = 299
    nc = 459
    for date in dates:
        for gkey in gkeys:
            outpath = (
                f'epa2022v1/{gkey}/{gkey}_{date:%Y-%m-%d}_epa2022v1_hc_22m.nc'
            )
            if os.path.exists(outpath):
                continue
            print(date, gkey)
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            try:
                gf = open_gdemis(date, gkey)
                to_ioapi(gf, outpath)
            except Exception as e:
                print(f'**WARNING:: Skipping {date} {gkey}: {e}')
                continue

        for pkey in pkeys:
            outpath = (
                f'epa2022v1/{pkey}/{pkey}_{date:%Y-%m-%d}_epa2022v1_hc_22m.nc'
            )
            if os.path.exists(outpath):
                continue
            print(date, pkey)
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            try:
                pf = open_ptemis(date, pkey)
                rpf = pt2gd(pf, nr, nc, ez=ez)
                to_ioapi(rpf, outpath)
            except IOError as e:
                print(f'**WARNING:: Skipping {date} {pkey}: {e}')
                continue
