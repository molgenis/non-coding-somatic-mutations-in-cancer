#!/usr/bin/bash

yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
}

## The file from which a function is called.
# source ${SCRIPT_PATH}variant_calling.sh