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
