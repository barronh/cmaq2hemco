__all__ = ['open_gdemis', 'open_ptemis']

import xarray as xr
xr.set_options(keep_attrs=True)


def open_gdemis(date, sector, cache=True):
    """
    Open gridded (gd) emission file for a date/sector combination

    Arguments
    ---------
    date : datetime
        Datetime object that supports strftime interface
    sector : str
        Named premerged sector
    cache : bool
        If True, store the input file locally, if false read in memory.

    Returns
    -------
    ef : xarray.Dataset
        xarray Dataset with lon, lat, and time meta variables.
    """
    from .utils import gd_file, open_date
    eroot = 'emis/2022v1/2022hc_cb6_22m/12US1/cmaq_cb6ae7'
    if sector == 'merged_nobeis_norwc':
        epath = (
            f'{eroot}/{sector}/emis_mole_all_%Y%m%d'
            + '_12US1_nobeis_norwc_2022hc_cb6_22m.ncf'
        )
    else:
        epath = (
            f'{eroot}/premerged_area/{sector}/emis_mole_{sector}_%Y%m%d'
            + '_12US1_cmaq_cb6ae7_2022hc_cb6_22m.ncf'
        )
    bucket = 'epa-2022-modeling-platform'
    try:
        ef = open_date(date, epath, bucket, cache=cache)
    except Exception as e:
        print(e)
        epath += '.gz'
        ef = open_date(date, epath, bucket, cache=cache)
    ef = gd_file(ef)

    return ef


def open_ptemis(date, sector, cache=True):
    """
    Open point (pt) emission file for a date/sector combination.
    This interface adds metadata from the stack groups file to
    create a single self-describing data file for subsequent processing.


    Arguments
    ---------
    date : datetime
        Datetime object that supports strftime interface
    sector : str
        Named point sector
    cache : bool
        If True, store the input file locally, if false read in memory.

    Returns
    -------
    ef : xarray.Dataset
        xarray Dataset with lon, lat, and time meta variables.
    """
    from .utils import se_file, open_date
    eroot = 'emis/2022v1/2022hc_cb6_22m/12US1/cmaq_cb6ae7'
    bucket = 'epa-2022-modeling-platform'
    try:
        spath = (
            f'{eroot}/{sector}/stack_groups_{sector}_%Y%m%d'
            + '_12US1_2022hc_cb6_22m.ncf'
        )
        sf = open_date(date, spath, bucket, cache=cache)
    except Exception:
        spath = (
            f'{eroot}/{sector}/stack_groups_{sector}_12US1_2022hc_cb6_22m.ncf'
        )
        sf = open_date(date, spath, bucket, cache=cache)

    epath = (
        f'{eroot}/{sector}/inln_mole_{sector}_%Y%m%d'
        + '_12US1_cmaq_cb6ae7_2022hc_cb6_22m.ncf'
    )
    ef = open_date(date, epath, bucket, cache=cache).isel(
        LAY=0, COL=0, drop=True
    )
    ef = se_file(sf, ef)
    return ef
