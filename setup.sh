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
    ./configure --prefix=${dir_packages}/python2
    make
    make test
    make install
    cd ..
    rm Python-2.7.15.tgz
    rm -rf Python-3.7.15
fi

# numpy

python2/bin/pip2.7 install numpy --no-deps

# astropy

python2/bin/pip2.7 install astropy --no-deps

# matplotlib

python2/bin/pip2.7 install matplotlib --no-deps

# astroML

python2/bin/pip2.7 install astroML


