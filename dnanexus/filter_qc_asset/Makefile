SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# samtools 1.3.1
	apt-get -q -y install liblzma-dev
	wget -nv https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
	tar -xvf samtools-1.3.1.tar.bz2
	cd samtools-1.3.1 && ./configure --prefix=/usr/local/ && make && sudo make install
	# bedtools 2.26
	wget -nv https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
	tar -xvf bedtools-2.26.0.tar.gz
	cd bedtools2 && make && sudo cp bin/* /usr/local/bin/
	# picard 1.92
	apt-get -q -y install ant openjdk-6-jdk
	wget -nv https://github.com/broadinstitute/picard/archive/1.92.zip
	unzip 1.92.zip
	cd picard-1.92 && ant -Djava6.home=/usr/lib/jvm/java-6-openjdk-amd64/ -Ddist=/picard -lib lib/ant package-commands
