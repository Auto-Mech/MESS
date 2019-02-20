# run with: . /path/to/fake-install.sh
export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export LD_LIBRARY_PATH=$THIS_DIR/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$THIS_DIR/lib:$LIBRARY_PATH
