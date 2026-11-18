# ============================================================================
# env.sh -- single source of paths for the SmbhBinary project.
#
#   source /path/to/smbh-binary/env.sh
#
# Machine-specific settings (scratch directory, SpaceHub location, SLURM account
# and partitions) go in env.local.sh next to this file.  It is gitignored; copy
# env.local.sh.example to start.  Anything already set in the environment wins.
#
# SLURM job scripts cannot locate themselves (sbatch runs a spooled copy), so the
# submit_*.sh scripts source this file and sbatch exports SMBH_ROOT to the job.
# ============================================================================

# locate this file (bash or zsh) -> SMBH_ROOT
_smbh_src=""
if [ -n "${BASH_SOURCE[0]:-}" ]; then _smbh_src="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then eval '_smbh_src="${(%):-%x}"'
fi
if [ -n "$_smbh_src" ] && [ -f "$(dirname "$_smbh_src")/env.sh" ]; then
    SMBH_ROOT="$(cd "$(dirname "$_smbh_src")" && pwd)"
fi
unset _smbh_src
if [ -z "${SMBH_ROOT:-}" ]; then
    echo "env.sh: cannot locate the repository; set SMBH_ROOT=/path/to/smbh-binary" >&2
    return 1 2>/dev/null || exit 1
fi
export SMBH_ROOT

# machine-specific overrides (gitignored)
[ -f "$SMBH_ROOT/env.local.sh" ] && . "$SMBH_ROOT/env.local.sh"

# SpaceHub headers (external library): default location inside the repo (gitignored)
if [ -z "${SPACEHUB:-}" ] && [ -f "$SMBH_ROOT/external/SpaceHub/src/spaceHub.hpp" ]; then
    SPACEHUB="$SMBH_ROOT/external/SpaceHub/src"
fi
export SPACEHUB

export SMBH_BIN="${SMBH_BIN:-$SMBH_ROOT/simulation/bin}"     # smbh_binary_PN / smbh_binary_Newt
export SMBH_DATA="${SMBH_DATA:-$SMBH_ROOT/data}"             # large data (gitignored)
export SMBH_SCRATCH="${SMBH_SCRATCH:-$SMBH_ROOT/scratch}"    # campaign run directories

# SLURM settings (optional; empty = cluster default)
export SMBH_SLURM_ACCOUNT="${SMBH_SLURM_ACCOUNT:-}"          # sbatch -A
export SMBH_SLURM_PARTITION="${SMBH_SLURM_PARTITION:-}"      # production jobs (~1 day, many cores)
export SMBH_SLURM_PARTITION_TEST="${SMBH_SLURM_PARTITION_TEST:-}"   # short test jobs

# slurm_submit <job.slurm> <log-dir> <partition> [extra sbatch args...]
#   Submits with log files in <log-dir> and the account from SMBH_SLURM_ACCOUNT.
slurm_submit() {
    local script="$1" logdir="$2" part="$3"; shift 3
    mkdir -p "$logdir"
    local args=(--export=ALL -o "$logdir/%x_%j.out" -e "$logdir/%x_%j.err")
    [ -n "$part" ] && args+=(-p "$part")
    [ -n "$SMBH_SLURM_ACCOUNT" ] && args+=(-A "$SMBH_SLURM_ACCOUNT")
    sbatch "${args[@]}" "$@" "$script"
}
