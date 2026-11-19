#!/usr/bin/env bash
# Shared environment loader for GWASfinemapping scripts (mirrors GWASprocessing/src/lib/env.sh).
# Safe to source from bash or zsh; detects cloud-synced folders and container runtimes.

if [[ -n "${BASH_VERSION:-}" ]]; then
  GWASFM_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
elif [[ -n "${ZSH_VERSION:-}" ]]; then
  GWASFM_LIB_DIR="$(cd "$(dirname "$0")" && pwd)"
else
  GWASFM_LIB_DIR="$(pwd)/scripts/lib"
fi
GWASFM_ROOT="$(cd "${GWASFM_LIB_DIR}/../.." && pwd)"
export GWASFM_ROOT

gwasfm_in_cloud_sync() {
  case "$1" in
    *OneDrive*|*CloudStorage*|*Dropbox*|*iCloud*|*"Google Drive"*) return 0 ;;
    *) return 1 ;;
  esac
}

gwasfm_in_container() {
  [[ "${GWASFM_IN_CONTAINER:-0}" == "1" || -f "/.dockerenv" ||
     -n "${SINGULARITY_NAME:-}" || -n "${APPTAINER_NAME:-}" ]]
}

gwasfm_load_dotenv() {
  local file="${1:-${GWASFM_ROOT}/.env}" line key value
  [[ -f "$file" ]] || return 0
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^[[:space:]]*(#|$) ]] && continue
    [[ "$line" == *=* ]] || continue
    key="${line%%=*}"
    key="${key#export }"
    key="${key//[[:space:]]/}"
    [[ "$key" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
    value="${line#*=}"
    value="${value%\"}"; value="${value#\"}"
    value="${value%\'}"; value="${value#\'}"
    if ! printenv "$key" >/dev/null 2>&1; then
      export "$key=$value"
    fi
  done < "$file"
}

gwasfm_resolve_paths() {
  gwasfm_load_dotenv "${GWASFM_ROOT}/.env"

  local default_cache="${HOME}/.cache/gwasfinemapping"
  if [[ -z "${GWASFM_WORK_ROOT:-}" ]]; then
    if gwasfm_in_cloud_sync "${GWASFM_ROOT}"; then
      export GWASFM_WORK_ROOT="${HOME}/GWASfinemapping-data/finemapping"
    else
      export GWASFM_WORK_ROOT="${GWASFM_ROOT}/finemapping"
    fi
  elif [[ "${GWASFM_WORK_ROOT}" != /* ]]; then
    export GWASFM_WORK_ROOT="${GWASFM_ROOT}/${GWASFM_WORK_ROOT}"
  fi

  if [[ -z "${GWASFM_TOOLS_ROOT:-}" ]]; then
    export GWASFM_TOOLS_ROOT="${default_cache}"
  elif [[ "${GWASFM_TOOLS_ROOT}" != /* ]]; then
    export GWASFM_TOOLS_ROOT="${GWASFM_ROOT}/${GWASFM_TOOLS_ROOT}"
  fi

  if [[ -z "${GWASFM_REFERENCE_ROOT:-}" ]]; then
    export GWASFM_REFERENCE_ROOT="${GWASFM_TOOLS_ROOT}/reference"
  fi
  if [[ -z "${GWASFM_CACHE_DIR:-}" ]]; then
    export GWASFM_CACHE_DIR="${GWASFM_TOOLS_ROOT}"
  fi
  if [[ -n "${GWASFM_TEMP_ROOT:-}" ]]; then
    if [[ "${GWASFM_TEMP_ROOT}" != /* ]]; then
      export GWASFM_TEMP_ROOT="${GWASFM_ROOT}/${GWASFM_TEMP_ROOT}"
    else
      export GWASFM_TEMP_ROOT
    fi
    export GWASFM_SCRATCH_DIR="${GWASFM_SCRATCH_DIR:-${GWASFM_TEMP_ROOT}/scratch}"
    export GWASFM_NXF_WORK="${GWASFM_NXF_WORK:-${GWASFM_TEMP_ROOT}/nextflow-work}"
    export GWASFM_TMP_DIR="${GWASFM_TMP_DIR:-${GWASFM_TEMP_ROOT}/tmp}"
  else
    export GWASFM_SCRATCH_DIR="${GWASFM_SCRATCH_DIR:-${GWASFM_CACHE_DIR}/scratch}"
    if [[ -z "${GWASFM_NXF_WORK:-}" ]]; then
      if gwasfm_in_cloud_sync "${GWASFM_ROOT}"; then
        export GWASFM_NXF_WORK="${GWASFM_TOOLS_ROOT}/nextflow-work"
      else
        export GWASFM_NXF_WORK="${GWASFM_TOOLS_ROOT}/nextflow-work"
      fi
    fi
  fi
  if [[ -z "${NXF_HOME:-}" ]] && gwasfm_in_cloud_sync "${GWASFM_ROOT}"; then
    export NXF_HOME="${GWASFM_TOOLS_ROOT}/nextflow-home"
  fi
  export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-${GWASFM_TOOLS_ROOT}/apptainer-cache}"
  export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-${APPTAINER_CACHEDIR}}"
}

gwasfm_activate_runtime() {
  gwasfm_resolve_paths
  export PYTHONNOUSERSITE=1
  export RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE
  export R_MAX_VSIZE="${R_MAX_VSIZE:-64Gb}"
  export OPENBLAS_NUM_THREADS="${GWASFM_BLAS_THREADS:-2}"
  export OMP_NUM_THREADS="${GWASFM_BLAS_THREADS:-2}"

  if gwasfm_in_container; then
    export GWASFM_IN_CONTAINER=1
    return 0
  fi

  local prefix="${GWASFM_TOOLS_ROOT}/env"
  if [[ -x "${prefix}/bin/python" ]]; then
    export GWASFM_ENV="${prefix}"
    case ":${PATH}:" in
      *":${prefix}/bin:"*) ;;
      *) export PATH="${prefix}/bin:${PATH}" ;;
    esac
    if [[ -d "${prefix}/lib/jvm" ]]; then
      export JAVA_HOME="${prefix}/lib/jvm"
    fi
  fi
  if [[ -d "${HOME}/.local/bin" ]] && ! command -v nextflow >/dev/null 2>&1; then
    export PATH="${HOME}/.local/bin:${PATH}"
  fi

  export RENV_PROJECT="${GWASFM_ROOT}"
  export RENV_PATHS_LIBRARY="${RENV_PATHS_LIBRARY:-${GWASFM_TOOLS_ROOT}/renv-library}"
  local _direct_r_lib="" _susie_dir=""
  if [[ -d "${RENV_PATHS_LIBRARY}" ]]; then
    _susie_dir="$(find "${RENV_PATHS_LIBRARY}" -mindepth 2 -maxdepth 6 \( -type d -o -type l \) -name "susieR" 2>/dev/null | head -n 1 || true)"
    [[ -n "${_susie_dir}" ]] && _direct_r_lib="$(dirname "${_susie_dir}")"
  fi
  if [[ -n "${_direct_r_lib}" && -e "${_direct_r_lib}/susieR/DESCRIPTION" ]]; then
    export R_LIBS_USER="${_direct_r_lib}"
    export R_PROFILE_USER=""
  elif [[ -f "${GWASFM_ROOT}/.Rprofile" ]]; then
    export R_PROFILE_USER="${GWASFM_ROOT}/.Rprofile"
  fi

  if [[ -z "${GWASFM_CONTAINER_IMAGE:-}" && -f "${GWASFM_TOOLS_ROOT}/gwasfinemapping.sif" ]]; then
    export GWASFM_CONTAINER_IMAGE="${GWASFM_TOOLS_ROOT}/gwasfinemapping.sif"
  fi
}

gwasfm_in_container_exec() {
  gwasfm_resolve_paths
  local engine="${GWASFM_CONTAINER:-apptainer}"
  local image="${GWASFM_CONTAINER_IMAGE:-${GWASFM_CONTAINER_REF:-ghcr.io/mattlee821/finemapping:latest}}"
  case "$engine" in
    docker)
      docker run --rm \
        -e GWASFM_IN_CONTAINER=1 \
        -e NXF_HOME=/tmp/gwasfm-nxf-home \
        -v "${GWASFM_ROOT}:${GWASFM_ROOT}" \
        -v "${GWASFM_REFERENCE_ROOT}:${GWASFM_REFERENCE_ROOT}" \
        -v "${GWASFM_CACHE_DIR}:${GWASFM_CACHE_DIR}" \
        -w "${GWASFM_ROOT}" \
        "$image" "$@"
      ;;
    apptainer|singularity)
      local bin
      bin="$(command -v apptainer || command -v singularity)"
      "$bin" exec \
        --env GWASFM_IN_CONTAINER=1 \
        --env NXF_HOME=/tmp/gwasfm-nxf-home \
        --bind "${GWASFM_ROOT}:${GWASFM_ROOT}" \
        --bind "${GWASFM_REFERENCE_ROOT}:${GWASFM_REFERENCE_ROOT}" \
        --bind "${GWASFM_CACHE_DIR}:${GWASFM_CACHE_DIR}" \
        --pwd "${GWASFM_ROOT}" \
        "$image" "$@"
      ;;
    *)
      echo "Unsupported GWASFM_CONTAINER engine: $engine (expected docker or apptainer)" >&2
      return 2
      ;;
  esac
}
