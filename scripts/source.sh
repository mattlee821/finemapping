#!/usr/bin/env bash
# Install the local tool environment and optionally pull the published tools image.
#
# Usage:
#   bash scripts/source.sh                    # tools and renv packages
#   bash scripts/source.sh --tools-only       # Hail, ldcov, Python, Java, and R
#   bash scripts/source.sh --r-packages       # restore the R project from renv.lock
#   bash scripts/source.sh --container apptainer
#   bash scripts/source.sh --container docker
#   bash scripts/source.sh --env-only
#   source scripts/source.sh                  # activate an installed environment

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${ROOT_DIR}/.env"
ENV_TEMPLATE="${ROOT_DIR}/.env.example"
CONTAINER_DEFAULT="ghcr.io/mattlee821/finemapping:latest"

load_env_defaults() {
  local file="$1" line key value
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

ensure_env_file() {
  if [[ ! -f "$ENV_FILE" ]]; then
    if [[ ! -f "$ENV_TEMPLATE" ]]; then
      echo "Missing environment template: $ENV_TEMPLATE" >&2
      return 1
    fi
    cp "$ENV_TEMPLATE" "$ENV_FILE"
    echo "Created .env from .env.example"
  fi
}

tools_root() {
  local root="${GWASFM_TOOLS_ROOT:-}"
  [[ -n "$root" ]] || root="${HOME}/.cache/gwasfinemapping"
  [[ "$root" == /* ]] || root="${ROOT_DIR}/${root}"
  printf '%s\n' "$root"
}

mamba_platform() {
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64) echo osx-arm64 ;;
    Darwin-x86_64) echo osx-64 ;;
    Linux-x86_64) echo linux-64 ;;
    Linux-aarch64|Linux-arm64) echo linux-aarch64 ;;
    *)
      echo "Unsupported platform $(uname -s)-$(uname -m)." >&2
      return 1
      ;;
  esac
}

install_micromamba() {
  local root bin platform
  root="$(tools_root)"
  if command -v micromamba >/dev/null 2>&1; then
    MAMBA_BIN="$(command -v micromamba)"
  else
    bin="${root}/bin/micromamba"
    platform="$(mamba_platform)"
    mkdir -p "$(dirname "$bin")"
    echo "Downloading micromamba for ${platform} into $(dirname "$bin")"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${platform}/latest" \
      | tar -xj -C "$root" bin/micromamba
    chmod +x "$bin"
    MAMBA_BIN="$bin"
  fi
  export MAMBA_ROOT_PREFIX="${root}/mamba"
}

file_hash() {
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$1" | awk '{print $1}'
  else
    sha256sum "$1" | awk '{print $1}'
  fi
}

install_mac_intel_lazybgen() {
  local python="$1"
  "$python" - <<'PY'
import hashlib
import json
import pathlib
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

version = "0.2.0"
expected_sha256 = "71f65d704b7326da55c45c63d9a956acb55b4b07fc07f5e6e4ac126ab916ab3e"
try:
    import lazybgen
    from importlib.metadata import version as package_version
    if package_version("lazybgen") == version:
        print(f"Using lazybgen {version}")
        raise SystemExit(0)
except ImportError:
    pass

with urllib.request.urlopen(f"https://pypi.org/pypi/lazybgen/{version}/json") as response:
    release = json.load(response)
sdist = next(item for item in release["urls"] if item["filename"].endswith(".tar.gz"))
if sdist["digests"]["sha256"] != expected_sha256:
    raise SystemExit("Unexpected lazybgen source archive checksum from PyPI.")

with tempfile.TemporaryDirectory(prefix="gwasfinemapping-lazybgen-") as temp:
    archive_path = pathlib.Path(temp) / sdist["filename"]
    with urllib.request.urlopen(sdist["url"]) as response:
        archive_bytes = response.read()
    if hashlib.sha256(archive_bytes).hexdigest() != expected_sha256:
        raise SystemExit("lazybgen source archive checksum verification failed.")
    archive_path.write_bytes(archive_bytes)
    with tarfile.open(archive_path, "r:gz") as archive:
        archive.extractall(temp, filter="data")

    source = pathlib.Path(temp) / f"lazybgen-{version}"
    setup = source / "setup.py"
    text = setup.read_text()
    marker = (
        '            if sys.platform == "darwin":\n'
        '                compile_args.extend(["-stdlib=libc++"])'
    )
    replacement = marker + (
        '\n                if platform.machine() == "x86_64":\n'
        '                    compile_args.append("-DZSTD_DISABLE_ASM")'
    )
    if text.count(marker) != 1:
        raise SystemExit("Could not uniquely locate lazybgen's macOS zstd build flags.")
    setup.write_text(text.replace(marker, replacement))
    subprocess.run(
        [sys.executable, "-m", "pip", "install", "--no-deps", str(source)],
        check=True,
    )
PY
}

install_python_tools() {
  local python="$1"
  if [[ "$(uname -s)-$(uname -m)" == "Darwin-x86_64" ]]; then
    install_mac_intel_lazybgen "$python"
  fi
  # PYTHONNOUSERSITE: pip must not treat packages in ~/.local as installed, because the environment is activated without them
  PYTHONNOUSERSITE=1 "$python" -m pip install --upgrade pip   # the pip that conda ships (24.x) cannot resolve the pinned stack
  PYTHONNOUSERSITE=1 "$python" -m pip install -r "${ROOT_DIR}/requirements-tools.txt"
}

install_environment() {
  local root prefix stamp spec requirements hash
  root="$(tools_root)"
  prefix="${root}/env"
  stamp="${prefix}/.tools-spec.sha256"
  spec="${ROOT_DIR}/environment.yml"
  requirements="${ROOT_DIR}/requirements-tools.txt"
  hash="$(file_hash "$spec").$(file_hash "$requirements")"
  install_micromamba

  if [[ "${FORCE:-0}" -ne 1 && -x "${prefix}/bin/python" &&
        -f "$stamp" && "$(cat "$stamp")" == "$hash" ]]; then
    echo "Tool environment is up to date: $prefix"
  elif [[ -d "${prefix}/conda-meta" ]]; then
    echo "Updating tool environment from environment.yml"
    "$MAMBA_BIN" install -y -q -p "$prefix" -f "$spec"
  else
    echo "Creating tool environment (first install may take several minutes)"
    "$MAMBA_BIN" create -y -q -p "$prefix" -f "$spec"
  fi
  install_python_tools "${prefix}/bin/python"
  mkdir -p "$prefix"
  printf '%s\n' "$hash" > "$stamp"
}

activate_tools() {
  local root prefix
  root="$(tools_root)"
  prefix="${root}/env"
  if [[ ! -x "${prefix}/bin/python" ]]; then
    echo "Tool environment not found. Run 'bash scripts/source.sh --tools-only' first." >&2
    return 1
  fi
  export GWASFM_TOOLS_ROOT="$root"
  export GWASFM_ENV="$prefix"
  export PATH="${prefix}/bin:${PATH}"
  export PYTHONNOUSERSITE=1
  if [[ -d "${prefix}/lib/jvm" ]]; then
    export JAVA_HOME="${prefix}/lib/jvm"
  fi
}

validate_tools() {
  python3 -c 'import hail, lazybgen, ldcov, pandas, pyarrow; from importlib.metadata import version; print("Hail", version("hail"), "| ldcov", version("ldcov"), "| lazybgen", version("lazybgen"))'
}

restore_r_packages() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is required to restore the R package environment." >&2
    return 1
  fi
  local library_root
  library_root="$(tools_root)/renv-library"
  mkdir -p "$library_root"
  (cd "$ROOT_DIR" && RENV_PATHS_LIBRARY="$library_root" RENV_CONFIG_PROMPT=FALSE \
    Rscript -e 'renv::restore(prompt = FALSE)')
}

pull_container() {
  local container="$1" root command sif
  root="$(tools_root)"
  case "$container" in
    docker)
      if ! command -v docker >/dev/null 2>&1; then
        echo "GWASFM_CONTAINER=docker but docker is not on PATH." >&2
        return 1
      fi
      docker pull "${GWASFM_CONTAINER_REF:-$CONTAINER_DEFAULT}"
      ;;
    apptainer|singularity)
      command="$(command -v apptainer || command -v singularity || true)"
      if [[ -z "$command" ]]; then
        echo "Apptainer/Singularity is not on PATH." >&2
        return 1
      fi
      sif="${root}/gwasfinemapping.sif"
      mkdir -p "$root"
      "$command" pull --force "$sif" "docker://${GWASFM_CONTAINER_REF:-$CONTAINER_DEFAULT}"
      export GWASFM_CONTAINER_IMAGE="$sif"
      ;;
    *)
      echo "Container must be docker or apptainer (got $container)." >&2
      return 2
      ;;
  esac
}

main() {
  set -euo pipefail
  local tools_only=0 r_only=0 env_only=0 container="" container_set=0
  FORCE=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --tools-only) tools_only=1 ;;
      --r-packages) r_only=1 ;;
      --env-only) env_only=1 ;;
      --force) FORCE=1 ;;
      --container)
        [[ $# -ge 2 ]] || { echo "--container requires docker or apptainer" >&2; return 2; }
        container="$2"
        container_set=1
        shift
        ;;
      --container=*) container="${1#*=}"; container_set=1 ;;
      -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}"; return 0 ;;
      *) echo "Unknown option: $1" >&2; return 2 ;;
    esac
    shift
  done

  ensure_env_file
  load_env_defaults "$ENV_FILE"
  if [[ "$container_set" -eq 0 ]]; then
    container="${GWASFM_CONTAINER:-}"
  fi
  if [[ "$env_only" -eq 1 ]]; then
    return 0
  fi
  if [[ "$r_only" -eq 1 ]]; then
    activate_tools
    restore_r_packages
  else
    install_environment
    activate_tools
    validate_tools
    [[ "$tools_only" -eq 1 ]] || restore_r_packages
  fi
  [[ -z "$container" ]] || pull_container "$container"
  echo "Setup complete. In a new shell, run 'source scripts/source.sh' to activate the tools."
}

if [[ "${BASH_SOURCE[0]}" != "$0" ]]; then
  load_env_defaults "$ENV_FILE"
  activate_tools
else
  main "$@"
fi
