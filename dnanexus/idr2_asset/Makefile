SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# Prerequisites for IDR on Ubuntu 14.04
	sudo apt-get -y install python3-dev=3.4.0-0ubuntu2
	sudo apt-get -y install python3-numpy=1:1.8.2-0ubuntu0.1
	sudo apt-get -y install python3-scipy=0.13.3-1build1
	sudo apt-get -y install python3-matplotlib=1.3.1-1ubuntu5
	sudo apt-get -y install python3-setuptools=3.3-1ubuntu2
	# IDR 2.0.4.1
	wget -nv https://github.com/kundajelab/idr/archive/2.0.4.1.tar.gz
	tar -xvf 2.0.4.1.tar.gz
	cd idr-2.0.4.1 && python3 setup.py install
