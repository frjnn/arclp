#!/bin/sh
# Downloads the Netlib .mat files that accompany the MATLAB arcLP code.
set -eu
dest="${1:-data/arcLP}"
if [ -d "$dest" ]; then
    echo "$dest already present"
    exit 0
fi
mkdir -p "$(dirname "$dest")"
git clone --depth 1 https://github.com/yaguangyang/arcLP.git "$dest"
