#! /bin/bash -f 

dir_packages="packages"


if [ ! -e $dir_packages ]
then
    mkdir $dir_packages
fi

cd $dir_packages

dir_packages=`pwd`

# python 3

if [ ! -e "python2" ]
then
    wget https://www.python.org/ftp/python/2.7.15/Python-2.7.15.tgz
    tar -xzf Python-2.7.15.tgz
    mkdir python2
    cd Python-2.7.15
    ./configure --prefix=${dir_packages}/python2 \
                LDFLAGS="-L${dir_packages}/python2/extlib/lib -Wl,--rpath=${dir_packages}/python2/lib -Wl,--rpath=${dir_packages}/python2/extlib/lib" \
                CPPFLAGS="-I${dir_packages}/python2/extlib/include"
    make
    make test
    make install
    cd ..
    rm Python-2.7.15.tgz
    rm -rf Python-2.7.15
fi

# install pip if it is not preinstalled.

python2/bin/python2.7 -m ensurepip --default-pip

# numpy

python2/bin/python2.7 -m pip install numpy

# astropy

python2/bin/python2.7 -m pip install astropy

# matplotlib

python2/bin/python2.7 -m pip installmatplotlib

# astroML

python2/bin/python2.7 -m pip install astroML_addons

# h5py

python2/bin/python2.7 -m pip install h5py

# emcee

python2/bin/pip2.7 python2.7 -m pip install emcee

# pandas

python2/bin/pip2.7 python2.7 -m pip install pandas

# make each directory

mkdir {data,output}
