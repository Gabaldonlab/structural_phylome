# install pdbx
git clone https://github.com/soedinglab/pdbx $CONDA_PREFIX/lib/pdbx
cd $CONDA_PREFIX/lib/pdbx
mkdir build
cd build
cmake -DUserInstallOption=ON ../
make install

cp -r build/lib/pdbx/ $CONDA_PREFIX/lib/python3.10/
