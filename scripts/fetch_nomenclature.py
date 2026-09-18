#!/usr/bin/env python3
"""Fetch pinned IPD-IMGT/HLA nomenclature from its authoritative repository."""
import hashlib
import json
from pathlib import Path
import urllib.request

COMMIT = '5b915f27f7f620361cf83cb626eeac8a03c0247c'
BASE = f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{COMMIT}/'
ROOT = Path(__file__).resolve().parents[1] / 'workflow'
FILES = {
    'wmda/hla_nom_g.txt': 'hla-analysis/source/imgt/wmda/hla_nom_g.txt',
    'Allelelist_history.txt': 'hla-bench/source/imgt/Allelelist_history.txt',
    'Allelelist.txt': 'hla-spechla-pg/source/Allelelist.txt',
    'LICENCE.md': 'hla-analysis/source/imgt/LICENCE.md',
}

def main():
    records = []
    for source, target in FILES.items():
        data = urllib.request.urlopen(BASE + source, timeout=60).read()
        path = ROOT / target
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        records.append(dict(path=target, url=BASE+source,
                            sha256=hashlib.sha256(data).hexdigest()))
    alias = ROOT / 'hla-spechla-pg/source/hla_nom_g.txt'
    alias.parent.mkdir(parents=True, exist_ok=True)
    alias.write_bytes((ROOT / FILES['wmda/hla_nom_g.txt']).read_bytes())
    records.append(dict(path=str(alias.relative_to(ROOT)), url=BASE+'wmda/hla_nom_g.txt',
                        sha256=hashlib.sha256(alias.read_bytes()).hexdigest()))
    (ROOT / 'nomenclature-downloads.json').write_text(json.dumps(records, indent=2)+'\n')
    print(f'Fetched {len(records)} files from IPD-IMGT/HLA v3.65.0 ({COMMIT[:8]}).')

if __name__ == '__main__':
    main()
