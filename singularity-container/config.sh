#!/bin/sh

srcdir=$(dirname "$0")
img=baseimage.sif
dstdir=${TMP:-/tmp}/pysecdec-github-runner

if [ ! -e "$srcdir/$img" ]; then
    echo "* building $img (requires root access)"
    $srcdir/build.sh || exit 1
fi
mkdir -p "$dstdir/work" "$dstdir/runner" "$dstdir/homedir"
if [ ! -e "$dstdir/$img" ]; then
    echo "* copying $img to $dstdir"
    cp -a "$srcdir/$img" "$dstdir/$img"
fi
cd "$dstdir"
if [ ! -e runner/config.sh ]; then
    echo "* fetching github runner"
    curl -o runner/runner.tgz -L https://github.com/actions/runner/releases/download/v2.277.1/actions-runner-linux-x64-2.277.1.tar.gz
    tar vxf runner/runner.tgz -C runner 
    rm runner/runner.tgz
    printf '#!/bin/sh\ncd /runner\nexec daemon -n runner -o /tmp/runner.log ./run.sh\n' >runner/service.sh
    chmod 755 runner/service.sh
fi
pyver=$(singularity exec -C $img python3 --version)
pyver=${pyver#Python }
pyver=python${pyver%.*}
singularity exec -C -W work -H homedir:/homedir -B runner:/runner --pwd /runner $img ./config.sh --replace true --work _work --labels $pyver,${img%.sif},ubuntu-latest
