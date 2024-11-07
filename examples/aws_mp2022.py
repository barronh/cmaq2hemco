import pandas as pd
import numpy as np
from cmaq2hemco.mp2022 import process_dates



debug = True
dates = pd.date_range('2022-01-01', '2022-12-31', freq='d')
# excluding burn and biogenics
# openburn
# beis4
gkeys = """rwc
rail
onroad_gas
onroad_diesel
onroad_ca_adj_gas
onroad_ca_adj_diesel
np_solvents
np_oilgas
nonroad_gas
nonroad_diesel
nonpt
mexico_onroad
livestock
canada_ptdust_adj
canada_onroad
canada_og2D
canada_afdust_adj
airports
afdust_adj""".split()

# ignoring fires to avoid duplication
# ptagfire
# ptfire-wild
# ptfire-rx
# ptfire_othna
pkeys = """ptegu
canmex_point
cmv_c1c2_12
cmv_c3_12
pt_oilgas
ptnonipm""".split()

elat = np.linspace(15, 65, 51)
elon = np.linspace(-130, -50, 81)
# allocating emissions to level based on stack height (STKHT) and Briggs
# plume rise using approximate level heights from 
# http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
# and assuming constant T=288.15 and P=101325
nz = 11
ez = np.array([
    -6,  123.,  254.,  387.,  521.,  657.,  795.,  934., 1075.,
    1218., 1363., 1510., 1659., 1860., 2118., 2382., 2654., 2932.,
    3219., 3665., 4132., 4623., 5142., 5692.
])[:nz]
if debug:
    dates = dates[:1]
    pkeys = pkeys[:1]
    gkeys = gkeys[:1]
    print('**WARNING: in debug mode; only processing')
    print(dates)
    print(pkeys)
    print(gkeys)

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

