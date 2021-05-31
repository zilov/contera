set -e

echo "Building contera..."

cp -r $SRC_DIR/* $PREFIX/

mkdir $PREFIX/bin
cd $PREFIX/bin
ln -s  $PREFIX/contera.py ./contera
chmod +x ./contera
