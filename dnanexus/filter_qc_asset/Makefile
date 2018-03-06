SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# needed for samtools
	apt-get -q -y install liblzma-dev
	# samtools 0.1.19
	# wget -nv https://github.com/samtools/samtools/archive/0.1.19.zip
	# unzip 0.1.19.zip
	# cd samtools-0.1.19 && make && cp samtools /usr/local/bin
	# samtools 1.2
	wget -nv https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
	tar -xvf samtools-1.2.tar.bz2
	cd samtools-1.2 && make && sudo make install
	# needed for picard
	apt-get -q -y install ant openjdk-6-jdk	
	# picard 1.92
	# wget -nv https://github.com/broadinstitute/picard/archive/1.92.zip
	# unzip 1.92.zip
	# cd picard-1.92 && ant -Djava6.home=/usr/lib/jvm/java-6-openjdk-amd64/ -Ddist=/picard -lib lib/ant package-commands
	# picard 1.126
	wget https://github.com/broadinstitute/picard/releases/download/1.126/picard-tools-1.126.zip
	unzip picard-tools-1.126.zip
	sudo cp -r picard-tools-1.126 /picard
