SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# needed for samtools
	apt-get install liblzma-dev
	# samtools 0.1.19
	# wget -nv https://github.com/samtools/samtools/archive/0.1.19.zip
	# unzip 0.1.19.zip
	# cd samtools-0.1.19 && make && cp samtools /usr/local/bin
	# samtools 1.2
	wget -nv https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
	tar -xvf samtools-1.2.tar.bz2
	cd samtools-1.2 && make && sudo make install
	# samtools 1.6
	# wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
	# tar -xvf samtools-1.6.tar.bz2
	# cd samtools-1.6 && ./configure --prefix=/usr/local/ && make && make install
