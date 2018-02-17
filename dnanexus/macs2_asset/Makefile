SHELL=/bin/bash -e

all:
	# Trust the signing key for this repo. Reference: http://cran.rstudio.com/bin/linux/ubuntu/README.html
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# macs2 2.1.1.20160309
	# sudo pip install -U --no-deps MACS2==2.1.1.20160309
	wget https://pypi.python.org/packages/9f/99/a8ac96b357f6b0a6f559fe0f5a81bcae12b98579551620ce07c5183aee2c/MACS2-2.1.1.20160309.tar.gz
	tar -xf MACS2-2.1.1.20160309.tar.gz
	cd MACS2-2.1.1.20160309 && sudo python setup.py install
