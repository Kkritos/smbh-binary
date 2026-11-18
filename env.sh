# ============================================================================
# env.sh -- single source of paths for the SmbhBinary project.
#
#   source /path/to/SmbhBinary/env.sh
#
# Works on the Mac and on SeaWulf.  SLURM job scripts cannot locate themselves
# (sbatch runs a spooled copy), so they read SMBH_ROOT from the environment
# exported by the submit_*.sh scripts, falling back to the SeaWulf clone path.
# Any variable already set in the environment wins over the defaults below.
# ============================================================================

if [ -n "${BASH_SOURCE[0]:-}" ] && [ -f "$(dirname "${BASH_SOURCE[0]}")/env.sh" ]; then
    SMBH_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
export SMBH_ROOT="${SMBH_ROOT:-/path/to/smbh-binary}"

# SpaceHub headers (external library, not part of this repo)
if [ -z "${SPACEHUB:-}" ]; then
    for d in "$HOME/codes/SpaceHub/src" "$HOME/codes/SpaceHub/src"; do
        [ -f "$d/spaceHub.hpp" ] && { SPACEHUB="$d"; break; }
    done
fi
export SPACEHUB

export SMBH_BIN="${SMBH_BIN:-$SMBH_ROOT/simulation/bin}"     # smbh_binary_PN / smbh_binary_Newt
export SMBH_DATA="${SMBH_DATA:-$SMBH_ROOT/data}"             # large data (gitignored)

# Scratch root for SLURM campaigns (SeaWulf).  Each campaign uses its own
# subdirectory, e.g. $SMBH_SCRATCH/smbh_bhstar_encounters.
if [ -z "${SMBH_SCRATCH:-}" ]; then
    if [ -d "/path/to/scratch/$USER" ]; then SMBH_SCRATCH="/path/to/scratch/$USER"
    else SMBH_SCRATCH="$SMBH_ROOT/scratch"; fi
fi
export SMBH_SCRATCH
