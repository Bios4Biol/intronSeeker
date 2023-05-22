#!/bin/bash --login
# The --login ensures the bash configuration is loaded,
# enabling Conda.

# Enable strict mode.
set -euo pipefail
# ... Run whatever commands ...

# Temporarily disable strict mode and activate conda:
set +euo pipefail
conda activate ISeeker_environment

# Re-enable strict mode:
set -euo pipefail

# exec the final command:
exec "$@"
