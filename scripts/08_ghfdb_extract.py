from config import *
"""Extract GHFDB (IHFC Release 2024, v2026.03) entries around Espírito Santo."""
import pandas as pd
from config import GHFDB_TXT
d = pd.read_csv(GHFDB_TXT, sep='\t', skiprows=12, encoding='latin-1', low_memory=False)
d = d.loc[:, ~d.columns.str.startswith('Unnamed')]
d['lat'] = pd.to_numeric(d.lat_NS, errors='coerce'); d['lon'] = pd.to_numeric(d.long_EW, errors='coerce')
d['q'] = pd.to_numeric(d.q, errors='coerce')
r = d[d.lat.between(-23, -15.5) & d.lon.between(-44, -38.5)]
r.to_csv('ghfdb_regiao.csv', index=False)
print(len(r), 'entries saved to ghfdb_regiao.csv')
