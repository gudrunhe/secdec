#!/bin/sh

# This script builds release files from the latest commit,
# ready to be uploaded to PyPI.
#
# To this end, it first clones the repository to a temporary
# location, creates a source distribution, and then builds this
# source distribution in multiple Singularity containers based
# on manylinux images.
#
# In the end the script runs checks on the built files.

c() {
    echo "# $@"
    "$@"
    if [ $? -ne 0 ]; then
        echo "# -> failed"
        exit 1
    fi
}

python=${PYTHON:-python3}
tmpdir=$(mktemp -d)

echo "### Installing Python dependencies"
c $python -m pip install auditwheel build twine

echo "### Cloning into $tmpdir"
c git clone . $tmpdir
echo "Latest commit:"
c git -C "$tmpdir" --no-pager log -1 --pretty=oneline

echo "### Building a source distribution in $tmpdir"
c $python -m build -s "$tmpdir"
name=$(ls "$tmpdir/dist/" | head -1)
name=${name%.tar.gz}
if [ ! -f "$tmpdir/dist/$name.tar.gz" ]; then
    echo "Can't find the distfile in $tmpdir/dist"
    exit 1
fi

c mkdir -p dist
c cp -ai "$tmpdir/dist/$name.tar.gz" dist/

for tag in manylinux2010_x86_64 manylinux2014_x86_64 manylinux_2_24_x86_64 manylinux_2_28_x86_64; do
    image="$tag.sif"
    image_src="docker://quay.io/pypa/$tag"

    if [ -e "dist/$name-$tag.whl" ]; then
        echo "File dist/$name-$tag.whl already exists; skip."
        continue
    fi

    if [ -e "$image" ]; then
        echo "### Using the existing $image (delete it to update it)"
    else
        echo "### Downloading $image"
        c singularity build "$image" "$image_src"
    fi

    echo "### Unpacking the distribution in $tmpdir/build-$tag"
    mkdir "$tmpdir/build-$tag"
    c tar -C "$tmpdir/build-$tag" -x -f "$tmpdir/dist/$name.tar.gz"

    dir="$tmpdir/build-$tag/$name"
    if [ ! -d "$dir" ]; then
        echo "Directory doesn't exist: $dir"
    fi

    echo "### Building $tag"
    c singularity run -C \
        -B "$dir":/src \
        -W "$TMP" \
        --pwd /src \
        "$image" \
        /opt/python/cp38-cp38/bin/python -m build .

    echo "### Saving the wheel from $dir/dist/ into dist/"
    c ls -l "$dir"/dist/
    bdist=$(ls "$dir/dist/$name-"*.whl | head -1)
    if [ ! -f "$bdist" ]; then
        echo "Can't find the distfile in $dir/dist"
        exit 1
    fi

    echo "The result is now at dist/$name-py3-none-$tag.whl"
    c cp -ai "$bdist" "dist/$name-py3-none-$tag.whl"

    echo "### Removing $dir"
    c rm -rf "$dir"
done

echo "### Removing $tmpdir"
rm c -rf "$tmpdir"

echo "### Double-check the release files by running:"
echo "$python -m twine check 'dist/$name.tar.gz'"

for w in "dist/$name"*.whl; do
    echo "$python -m twine check '$w'"
    echo "$python -m auditwheel show '$w'"
done

echo "### If all the checks pass, then upload to PyPI by running:"
echo "$python -m twine upload 'dist/$name.tar.gz'"
for w in "dist/$name"*.whl; do
    echo "$python -m twine upload '$w'"
done
