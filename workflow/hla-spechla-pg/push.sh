#!/bin/bash
# Upload all pipeline scripts and job files to the cluster (streamed tar over the two-hop ssh transport).
set -euo pipefail
cd "$(dirname "$0")"
tar czf - --transform 's#^#x/#' *.py *.sh jobs | ssh -o BatchMode=yes ddbj "ssh -o BatchMode=yes -o IdentitiesOnly=yes -i ~/.ssh/id_ed25519 a001 'cd /home/leechuck/hla/spechla-pg && rm -rf x && tar xzf - && mkdir -p scripts jobs && mv x/*.py x/*.sh scripts/ && cp x/jobs/* jobs/ && rm -rf x && chmod +x scripts/*.sh && ls scripts jobs'"
