# run with: . /path/to/fake-install.sh
export THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
export PROJECT_ROOT=$THIS_DIR/..

export PATH=$THIS_DIR/bin:$PATH
