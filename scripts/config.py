"""Paths to the input data. Edit DATA, EMAG2_TIF and GHFDB_TXT for your machine.

DATA must point to a folder with:
  GM_K_Percetual.grd, GM_TH_PPM.grd, GM_U_PPM.grd (+ .gi)   Project 1093 gamma-ray grids (50 m)
  M_SemIGRF.grd (+ .gi)                                      Project 1093 magnetic anomaly, IGRF removed
  Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp
  Espirito Santo/Estrutural/EspiritoSanto_estrutura.shp
  Espirito Santo/Litologia/EspiritoSanto_lito.shp
"""
import os
DATA = os.environ.get('ES_DATA', '../data/') .rstrip('/') + '/'
EMAG2_TIF = os.environ.get('EMAG2_TIF', DATA + 'EMAG2_V3_20170530_UpCont.tif')
GHFDB_TXT = os.environ.get('GHFDB_TXT', DATA + 'IHFC_2024_GHFDB_v.2026.03.txt')
