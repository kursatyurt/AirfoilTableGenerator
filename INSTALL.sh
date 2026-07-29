#!/usr/bin/env bash
# Sets up everything under ./opt and ./.venv. No sudo. Re-runnable.
set -euo pipefail
cd "$(dirname "$0")"
ROOT=$PWD
NP=$(getconf _NPROCESSORS_ONLN)

# ---- preflight ------------------------------------------------------------
missing=()
for c in git c++ mpicxx mpirun python3; do
  command -v "$c" >/dev/null || missing+=("$c")
done
if [ ${#missing[@]} -ne 0 ]; then
  echo "Missing: ${missing[*]}" >&2
  echo "Debian/Ubuntu: sudo apt install -y git build-essential libopenmpi-dev openmpi-bin python3 python3-venv" >&2
  echo "Fedora/RHEL:   sudo dnf install -y git gcc-c++ openmpi-devel python3" >&2
  exit 1
fi
CPU_ARCH=${SU2_CPU_ARCH:-native}
SU2_VERSION=${SU2_VERSION:-v8.5.0}
MACHINE_ID="$(uname -m)|$(sysctl -n hw.model 2>/dev/null || uname -n)|$(mpicxx --version | head -1)"
BUILD_STAMP=$(printf 'su2=%s\narch=%s\nmachine=%s' "$SU2_VERSION" "$CPU_ARCH" "$MACHINE_ID")
python3 -c 'import sys; sys.exit(sys.version_info < (3, 8))' \
  || { echo "need python >= 3.8, got $(python3 -V)" >&2; exit 1; }
free=$(df -Pk . | awk 'NR==2{print int($4/1048576)}')
[ "$free" -ge 8 ] || echo "WARNING: only ${free}G free here; the SU2 build wants ~8G." >&2

# ---- python env -----------------------------------------------------------
echo "==> python venv"
if [ ! -x .venv/bin/python ]; then
  # ponytail: stock Ubuntu python3 ships without ensurepip, so venv fails there.
  # virtualenv --user is the no-sudo way out; apt install python3-venv also works.
  python3 -m venv .venv >/dev/null 2>&1 \
    || { echo "    (python3-venv unavailable, falling back to virtualenv)";
         python3 -m virtualenv .venv 2>/dev/null \
         || { pip3 install --user -q virtualenv && python3 -m virtualenv .venv; }; }
fi
./.venv/bin/python -m pip install -q --upgrade pip
./.venv/bin/python -m pip install -q gmsh numpy scipy pandas matplotlib

mkdir -p opt

# ---- FALCON (mesher + airfoil database) -----------------------------------
echo "==> FALCON"
[ -d opt/FALCON/.git ] || git clone --depth 1 https://github.com/Prisha22/FALCON opt/FALCON
# ponytail: only meshing.py / read_airfoil.py / Airfoil_DAT_Selig are used. Its
# requirements.txt pulls the Tk GUI stack, so it is deliberately not installed.

# ---- SU2 ------------------------------------------------------------------
if [ -x opt/su2/bin/SU2_CFD ] && [ -f opt/su2/.build-machine ] && \
   [ "$(cat opt/su2/.build-machine)" = "$BUILD_STAMP" ]; then
  echo "==> SU2 release build matches this machine ($CPU_ARCH)"
else
  echo "==> SU2 native release build with MPI (~20-40 min on $NP cores)"
  if [ -d opt/SU2-src/.git ]; then
    git -C opt/SU2-src fetch --depth 1 origin "refs/tags/$SU2_VERSION:refs/tags/$SU2_VERSION"
    git -C opt/SU2-src checkout --detach "$SU2_VERSION"
  else
    git clone --depth 1 --branch "$SU2_VERSION" https://github.com/su2code/SU2 opt/SU2-src
  fi
  cd opt/SU2-src
  # SU2 bundles meson and ninja as submodules; meson.py initialises them.
  ./meson.py build -Dwith-mpi=enabled -Denable-autodiff=false \
                   -Denable-pywrapper=false -Dcpu-arch="$CPU_ARCH" \
                   --buildtype=release --prefix="$ROOT/opt/su2"
  NINJA=./ninja; [ -x "$NINJA" ] || NINJA="$ROOT/.venv/bin/ninja"
  [ -x "$NINJA" ] || { "$ROOT/.venv/bin/python" -m pip install -q ninja; NINJA="$ROOT/.venv/bin/ninja"; }
  "$NINJA" -C build install -j "$NP"
  cd "$ROOT"
  printf '%s\n' "$BUILD_STAMP" > opt/su2/.build-machine
fi

cat > env.sh <<EOF
export SU2_HOME=$ROOT/opt/SU2-src
export SU2_RUN=$ROOT/opt/su2/bin
export PATH=\$SU2_RUN:$ROOT/.venv/bin:\$PATH
export PYTHONPATH=$ROOT/opt/FALCON:\${PYTHONPATH:-}
EOF

# ---- verify ---------------------------------------------------------------
echo "==> verifying"
# shellcheck disable=SC1091
source ./env.sh
command -v SU2_CFD >/dev/null || { echo "SU2_CFD not on PATH after install" >&2; exit 1; }
python -c "import gmsh, numpy, scipy, pandas, matplotlib, read_airfoil" \
  || { echo "python imports failed" >&2; exit 1; }
python "$ROOT/mesh.py"
python "$ROOT/polar.py" --selftest

echo
echo "OK. Next:"
echo "  source env.sh"
echo "  AIRFOIL=naca0012 MACH=\"0.15\" RE=1e6 AOA=-4:16:2 TRANSITION=none bash tools/run_rotor_table.sh"
