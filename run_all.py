"""Run the whole es-crustal-thermal-state workflow (scripts 01-13) in order.

Usage (from the repository folder):
    python run_all.py              # run everything
    python run_all.py --from 9     # resume from script 09
    python run_all.py --only 12    # run a single script

Edit the three paths below for your machine, or set the environment
variables ES_DATA, EMAG2_TIF and GHFDB_TXT before running.
If GHFDB_TXT does not exist, the GHFDB .txt is extracted from the zip
given in GHFDB_ZIP into data/ghfdb/.
"""
import argparse
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path

# ----------------------------------------------------------------- paths
ES_DATA = os.environ.get('ES_DATA', r'C:/Mestrado/ES_DADOS')
EMAG2_TIF = os.environ.get(
    'EMAG2_TIF',
    r'C:/Mestrado/ES_DADOS/Dados Magnetométricos/EMAG2_V3_20170530_UpCont.tif')
GHFDB_ZIP = os.environ.get(
    'GHFDB_ZIP', r'C:/Mestrado/ES_DADOS/GHFDB/GHFBD-R2024_v.2026-03.zip')
GHFDB_TXT = os.environ.get(
    'GHFDB_TXT', r'C:/Mestrado/ES_DADOS/GHFDB/IHFC_2024_GHFDB_v.2026.03.txt')
# -----------------------------------------------------------------------

ROOT = Path(__file__).resolve().parent
ORDER = [
    '01_heat_production_50m.py',
    '02_colatina_corridor.py',
    '03_heat_production_125m_and_reference_geotherm.py',
    '04_heat_production_by_lithology.py',
    '05_emag2v3_crop.py',
    '06_aeromagnetic_fractal_fit.py',
    '07_centroid_estimates.py',
    '08_ghfdb_extract.py',
    '09_monte_carlo.py',
    '10_sensitivity.py',
    '11_figures_1to3.py',
    '12_figure_4.py',
    '13_figure_5.py',
]
SCRIPTS = ROOT / 'scripts'
(ROOT / 'results').mkdir(exist_ok=True)
(ROOT / 'figures').mkdir(exist_ok=True)


def ghfdb_txt():
    """Return the path to the GHFDB .txt, extracting it from the zip if needed."""
    if GHFDB_TXT and Path(GHFDB_TXT).is_file():
        return Path(GHFDB_TXT)
    out = ROOT / 'data' / 'ghfdb'
    found = sorted(out.glob('*.txt')) if out.is_dir() else []
    if found:
        return found[0]
    z = Path(GHFDB_ZIP)
    if not z.is_file():
        sys.exit(f'GHFDB not found: set GHFDB_TXT or GHFDB_ZIP (tried {z})')
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(z) as zf:
        txt = [n for n in zf.namelist() if n.lower().endswith('.txt')]
        if not txt:
            sys.exit(f'No .txt inside {z}')
        name = max(txt, key=lambda n: zf.getinfo(n).file_size)  # the data table
        target = out / Path(name).name
        with zf.open(name) as src, open(target, 'wb') as dst:
            dst.write(src.read())
    print(f'GHFDB extracted to {target}')
    return target


def check_inputs(ghfdb):
    data = Path(ES_DATA)
    need = [
        data / 'GM_K_Percetual.grd', data / 'GM_TH_PPM.grd', data / 'GM_U_PPM.grd',
        data / 'M_SemIGRF.grd',
        data / 'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp',
        data / 'Espirito Santo/Estrutural/EspiritoSanto_estrutura.shp',
        data / 'Espirito Santo/Litologia/EspiritoSanto_lito.shp',
        Path(EMAG2_TIF), ghfdb,
    ]
    missing = [str(p) for p in need if not p.is_file()]
    if missing:
        print('Missing input files:')
        for m in missing:
            print('  ', m)
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--from', dest='start', type=int, default=1)
    ap.add_argument('--only', type=int)
    args = ap.parse_args()

    if sys.version_info < (3, 10):
        sys.exit('Python 3.10 or newer is required (3.13 was used for the paper).')

    ghfdb = ghfdb_txt()
    check_inputs(ghfdb)

    env = os.environ.copy()
    env['ES_DATA'] = ES_DATA.replace('\\', '/')
    env['EMAG2_TIF'] = EMAG2_TIF.replace('\\', '/')
    env['GHFDB_TXT'] = str(ghfdb).replace('\\', '/')
    env['MPLBACKEND'] = 'Agg'          # no windows pop up for the figures
    env['PYTHONIOENCODING'] = 'utf-8'  # Portuguese messages print fine

    scripts = [SCRIPTS / n for n in ORDER]
    if args.only:
        scripts = [s for s in scripts if int(s.name[:2]) == args.only]
    else:
        scripts = [s for s in scripts if int(s.name[:2]) >= args.start]

    log = open(ROOT / 'run_all.log', 'w', encoding='utf-8')
    t_all = time.time()
    for s in scripts:
        print(f'\n=== {s.name} ===', flush=True)
        log.write(f'\n=== {s.name} ===\n')
        t0 = time.time()
        p = subprocess.run([sys.executable, str(s)], cwd=ROOT / 'results', env=env,
                           capture_output=True, text=True, encoding='utf-8',
                           errors='replace')
        print(p.stdout, end='')
        log.write(p.stdout)
        if p.returncode != 0:
            print(p.stderr)
            log.write(p.stderr)
            log.close()
            print(f'\nFAILED at {s.name}. Full log: run_all.log')
            print(f'After fixing, resume with: python run_all.py --from {int(s.name[:2])}')
            sys.exit(p.returncode)
        msg = f'--- {s.name} done in {time.time() - t0:.0f} s'
        print(msg)
        log.write(msg + '\n')

    msg = f'\nAll {len(scripts)} scripts finished in {(time.time() - t_all) / 60:.1f} min.'
    print(msg)
    print('Outputs: results/ (numbers, GeoTIFFs) and figures/ (paper figures).')
    print('Compare with the paper: A0 median ~1.99 uW/m3, q0 ~30 mW/m2, T(Moho) ~280 C.')
    log.write(msg + '\n')
    log.close()


if __name__ == '__main__':
    main()
