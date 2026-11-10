#!/usr/bin/env bash
# Copy the barcode sets shipped in the repository (frontend/bc/) into the app's barcode
# volume, once, after the first `docker compose up`. The app lists barcode sets from that
# volume (TBRC_BC_DIR), which starts empty. Skips if the volume already has content, so it
# never overwrites sets added through the app.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../.." && pwd)"
cd "${HERE}"

existing="$(docker compose exec -T tbrc sh -c "ls -A /var/lib/tbrc/bc 2>/dev/null | wc -l" | tr -d '[:space:]')"
if [[ "${existing}" != "0" ]]; then
    echo "skip: /var/lib/tbrc/bc already has ${existing} entries"
    exit 0
fi
docker compose cp "${REPO_ROOT}/frontend/bc/." "tbrc:/var/lib/tbrc/bc/"
echo "seeded barcode sets: $(docker compose exec -T tbrc sh -c 'ls /var/lib/tbrc/bc' | tr '\n' ' ')"
