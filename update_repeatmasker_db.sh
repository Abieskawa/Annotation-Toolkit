#!/usr/bin/env bash
set -euo pipefail

OLD_DB=$HOME/output/RepeatMasker_DB_Dfam38_RepBase20181026_BACKUP
NEW_DB=$HOME/output/RepeatMasker_DB
IMAGE=dfam/tetools:1.99
BASE=https://www.dfam.org/releases/Dfam_3.9/families/FamDB

rm -rf "$NEW_DB"
mkdir -p "$NEW_DB"

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$NEW_DB:/host_libraries" \
  "$IMAGE" \
  bash -lc 'cp -a /opt/RepeatMasker/Libraries/. /host_libraries/'

mkdir -p "$NEW_DB/famdb"
cd "$NEW_DB/famdb"
rm -f dfam*.h5 dfam*.h5.gz dfam*.md5 rmlib.config

for p in 0 4 10 11 12
do
  wget -q -c "$BASE/dfam39_full.${p}.h5.gz"
  wget -q -c "$BASE/dfam39_full.${p}.h5.gz.md5"
done

md5sum -c *.md5
gunzip -f dfam39_full.*.h5.gz

docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "$NEW_DB:/opt/RepeatMasker/Libraries" \
  -v "$OLD_DB:/old_db:ro" \
  "$IMAGE" \
  bash -lc '
    set -euo pipefail
    cd /opt/RepeatMasker
    tar -xzf /old_db/RepBaseRepeatMaskerEdition-20181026.tar.gz
    ./tetoolsDfamUpdate.pl
  '