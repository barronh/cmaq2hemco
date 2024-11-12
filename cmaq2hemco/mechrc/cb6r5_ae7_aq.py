__all__ = ['writeconfig']


cq2gc = {
    'ACET': [['ACET', '26/213/1007']],
    'ALD2': [['ALD2', '26/213/1007']],
    'ALDX': [['RCHO', '26/213/1007']],
    'ACROLEIN': [['MACR', '26/213/1007']],
    'BENZ': [['BENZ', '26/213/1007']],
    'ECH4': [['CH4', '1007']],
    'CL': [['Cl', '1007']],
    'CL2': [['Cl2', '1007']],
    'CLO': [['ClO', '1007']],
    'CLNO2': [['ClNO2', '1007']],
    'CLNO3': [['ClNO3', '1007']],
    'CO': [['CO', '26/211/1007']],
    'ETOH': [['EOH', '26/213/1007']],
    'ETHA': [['C2H6', '26/217/1007']],
    'ETH': [['C2H4', '26/213/1007']],
    'FACD': [['HCOOH', '26/213/1007']],
    'FORM': [['CH2O', '26/213/1007']],
    'GLY': [['GLYX', '26/214/1007']],  # Following MEK
    'GLYD': [['GLYC', '26/214/1007']],  # Following MEK
    'HCL': [['HCl', '1007']],
    'HOCL': [['HOCl', '1007']],
    'HONO': [['HNO2', '25/210/1007']],
    'HPLD': [['HPALD', '26/214/1007']],  # Following MEK
    'IOLE': [['PRPE', '26/215/1007']],
    'KET': [['MEK', '26/214/1007']],
    'MEPX': [['MP', '26/213/1007']],  # Following MEOH
    'MEOH': [['MOH', '26/213/1007']],
    'NH3': [['NH3', '26/213/1007']],
    'NO': [['NO', '115/25/210/1007']],
    'NO2': [['NO2', '25/210/1007']],
    'OLE': [['PRPE', '26/215/1007']],
    'PACD': [['MAP', '26/213/1007']],  # Following MEOH
    'PAR': [['ALK4', '26/212/1007']],
    'PNA': [['HNO4', '26/213/1007']],  # Following MEOH
    'PRPA': [['C3H8', '26/216/1007']],
    'SO2': [['SO2', '26/218/1007']],
    'SULF': [['SO4', '26/218/1007']],
    'TOL': [['TOLU', '26/213/1007']],
    'PEC': [['BCPI', '26/221/1007/70'], ['BCPO', '26/221/256/1007/71']],
    'PFE': [['pFe', '26/219/1007']],
    'POC': [['OCPI', '26/222/1007/72'], ['OCPO', '26/222/257/1007/73']],
    'PNH4': [['NH4', '26/218/1007']],  # need to confirm this
    'PNO3': [['NIT', '26/218/1007']],  # need to confirm this
    'PSO4': [['SO4', '26/219/1007']],
    'XYLMN': [['XYLE', '26/213/1007']],
}
# Ignore special species
for key in ['ALD2_PRIMARY', 'ALD2_PRIMARY', 'NH3_FERT']:
    cq2gc[key] = []
# ignore inventory meta variables
for key in ['HFLUX', 'VOC_INV', 'VOC_BEIS', 'CO2_INV', 'N2O_INV', 'CH4_INV']:
    cq2gc[key] = []
# Ignore duplicative species
# TOLU in TOL
# NOX = NO + NO2 + HONO
# CH4 in ECH4
for key in ['CH4', 'TOLU', 'NOX']:
    cq2gc[key] = []

# Ignore species GC doesn't have
for key in ['BUTADIENE13', 'ETHY', 'NAPH', 'UNR', 'UNK', 'NR', 'NVOL']:
    cq2gc[key] = []

# Ignore unused particulate species
for key in ['SOAALK', 'PAL', 'PCA', 'PCL', 'PFE', 'PH2O', 'PK', 'PSI', 'PTI']:
    cq2gc[key] = []

for key in ['PM2_5', 'PMC', 'PMG', 'PMN', 'PMOTHR', 'PNCOM']:
    cq2gc[key] = []


def writeconfig(outpath, year, sector, filepatt):
    from .core import writeconfig as genwriteconfig
    return genwriteconfig(outpath, year, sector, filepatt, cq2gc)
