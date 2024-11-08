import pandas as pd
import numpy as np
import os
from cmaq2hemco.mp2022 import open_gdemis, open_ptemis
from cmaq2hemco.utils import pt2hemco, gd2hemco, gd2matrix, gd2hemco_fast


debug = True
dates = pd.date_range('2022-01-01', '2022-12-31', freq='d')
# Define grid by edges
elat = np.linspace(15, 65, 501)
elon = np.linspace(-135, -50, 851)
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
gkeys = ['merged_nobeis_norwc']
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

if debug:
    dates = dates[:1]
    pkeys = pkeys[:0]
    gkeys = gkeys[:1]
    print('**WARNING: in debug mode; only processing')
    print(dates)
    print(pkeys)
    print(gkeys)

# if not os.path.exists('epa2022v1/matrix.csv'):
#     gf = open_gdemis(dates[0], gkeys[0])
#     mdf = gd2matrix(gf, elat, elon)
#     mdf.drop('geometry', axis=1).to_csv('epa2022v1/matrix.csv')
# matrix = pd.read_csv('epa2022v1/matrix.csv', index_col=['ROW', 'COL', 'lati', 'loni'])
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
            # open_gdemis could fail if date is not available (eg., only weekday)
            gf = open_gdemis(date, gkey)
        except Exception as e:
            print(f'**WARNING:: Skipping {date} {gkey}: {e}')
            continue
        # using bilinear interpolation of fluxes
        rgf = gd2hemco_fast(outpath, gf, elat, elon)
        # use matrix interpolation for fractional area overlap (slow)
        # rgf = gd2hemco(outpath, gf, elat, elon, matrix=matrix)
        del rgf, gf
    
    for pkey in pkeys:
        outpath = (
            f'epa2022v1/{pkey}/{pkey}_{date:%Y-%m-%d}_epa2022v1_hc_22m.nc'
        )
        if os.path.exists(outpath):
            continue
        print(date, pkey)
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        try:
            # open_ptemis could fail if date is not available (eg., only weekday)
            pf = open_ptemis(date, pkey)
        except IOError as e:
            print(f'**WARNING:: Skipping {date} {pkey}: {e}')
            continue
        rpf = pt2hemco(outpath, pf, elat, elon)  # apply plume rise
        del rpf, pf

