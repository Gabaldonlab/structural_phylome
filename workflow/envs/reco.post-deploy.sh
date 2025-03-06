# install Notung
wget https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip -O $CONDA_PREFIX/lib/Notung-2.9.1.5.zip
unzip $CONDA_PREFIX/lib/Notung-2.9.1.5.zip  -d $CONDA_PREFIX/lib/

echo -e '#!/usr/bin/env bash\njava -jar $CONDA_PREFIX/lib/Notung-2.9.1.5.jar $@' > $CONDA_PREFIX/bin/notung
chmod +x $CONDA_PREFIX/bin/notung


git clone https://github.com/JSdoubleL/DISCO $CONDA_PREFIX/lib/DISCO
sed -i '1s/^/\#!\/usr\/bin\/env python3\n/' $CONDA_PREFIX/lib/DISCO/disco.py

chmod +x $CONDA_PREFIX/lib/DISCO/disco.py
cp $CONDA_PREFIX/lib/DISCO/disco.py $CONDA_PREFIX/bin


# install astral-pro
git clone https://github.com/chaoszhang/ASTER $CONDA_PREFIX/lib/ASTER
cd $CONDA_PREFIX/lib/ASTER
make astral-pro
cp bin/astral-pro $CONDA_PREFIX/bin/
cd -
