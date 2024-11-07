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

process_dates(dates, gkeys, pkeys, elat, elon, ez, debug=True)
