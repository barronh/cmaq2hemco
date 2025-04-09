__all__ = ['writeconfig']


def writeconfig(outpath, year, sector, filepatt, cq2gc):
    import pandas as pd
    import xarray as xr
    from warnings import warn
    import os
    defaults = set()
    ignores = set()
    exdates = pd.date_range(f'{year}-01-01', freq='1d', periods=365)
    hcpatt = filepatt.replace('%Y', '$YYYY')
    hcpatt = hcpatt.replace('%m', '$MM')
    hcpatt = hcpatt.replace('%d', '$DD')
    for date in exdates:
        hcpath = date.strftime(filepatt)
        if os.path.exists(hcpath):
            break
    else:
        raise IOError('file not found')

    with open(outpath, 'w') as hcf:
        print(sector, hcpatt, end='', flush=True)
        hcfile = xr.open_dataset(hcpath)
        for cqkey, v in hcfile.data_vars.items():
            if cqkey in hcfile.sizes or cqkey in ('hyai', 'hybi'):
                continue
            elif cqkey in ('TOLU',):
                warn('TOLU mass is duplicated by TOL')
            if cqkey in cq2gc:
                gctrans = cq2gc.get(cqkey)
                if len(gctrans) == 0:
                    ignores.add(cqkey)
            else:
                defaults.add(cqkey)
                gctrans = [[cqkey, '1007']]
            for gckey, scale in gctrans:
                # if gckey in [
                #     'ACET', 'MEK', 'ALD2', 'PRPE', 'PRPA', 'BENZ', 'TOLU',
                #     'XYLE', 'EOH', 'ALK4', 'ISOP'
                # ]:
                #     units = 'kgC/m2/s'
                # else:
                # All variables are being entered as kg/m2/s, which HEMCO
                # should read correctly.
                units = v.units.strip()
                opts = dict(
                    unit=units,
                    gckey=gckey,
                    cqkey=cqkey,
                    sector=sector,
                    path=hcpatt,
                    scale=scale,
                    cat='1/2',
                    hier=50
                )
                hcf.write((
                    '0 EPA{year}_{gckey}__{sector}{cqkey} {path}  {cqkey}'
                    '       {year}-{year}/1-12/1-31/0-24 C xyz  {unit}'
                    '  {gckey}   {scale}     {cat} {hier}\n'
                ).format(year=year, **opts))
        print()
    print('Ignored', sorted(ignores))
    print('Defaults', sorted(defaults))
