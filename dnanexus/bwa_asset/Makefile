SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# bwa
	wget https://github.com/lh3/bwa/releases/download/v0.7.13/bwakit-0.7.13_x64-linux.tar.bz2
	bzip2 -d bwakit-0.7.13_x64-linux.tar.bz2
	tar -xvf bwakit-0.7.13_x64-linux.tar
	cd bwa.kit && sudo cp bwa /usr/local/bin
