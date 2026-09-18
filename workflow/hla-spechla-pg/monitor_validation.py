#!/usr/bin/env python3
"""Collect and score the fixed cohort until completion or the prespecified deadline."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from fetch_validation import fetch
from score_validation import evaluate
from score_comparators import evaluate as evaluate_comparators


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('release'); p.add_argument('--out',type=Path,required=True)
    p.add_argument('--deadline',default='2026-09-18T09:02:07+00:00')
    p.add_argument('--interval',type=int,default=180)
    a=p.parse_args(); deadline=datetime.fromisoformat(a.deadline)
    if deadline.tzinfo is None or a.interval<30: p.error('Timezone and interval >=30 seconds required')
    cohort=Path(__file__).parent/'validation/20260918/cohort.tsv'
    a.out.mkdir(parents=True,exist_ok=True)
    while True:
        now=datetime.now(timezone.utc)
        at_deadline=now>=deadline
        try:
            fetch(a.release,a.out/'snapshot')
            evaluate(a.out/'snapshot',cohort,a.out/'analysis')
            evaluate_comparators(a.out/'snapshot',cohort,a.out/'analysis')
            record=json.loads((a.out/'analysis/analysis.json').read_text())
            record.update(observed_at=now.isoformat(),deadline=deadline.isoformat(),deadline_reached=at_deadline)
            (a.out/'monitor.json').write_text(json.dumps(record,indent=2)+'\n')
            if record['complete'] or at_deadline:
                return
        except Exception as exc:
            with open(a.out/'monitor-errors.log','a') as f:
                f.write(f'{now.isoformat()} {exc!r}\n')
            if at_deadline:
                raise
        time.sleep(min(a.interval,max(1,(deadline-datetime.now(timezone.utc)).total_seconds())))


if __name__=='__main__': main()
