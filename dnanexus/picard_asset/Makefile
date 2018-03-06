SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# needed for picard
	sudo apt-get -y install ant openjdk-6-jdk	
	# picard 1.92
	# wget -nv https://github.com/broadinstitute/picard/archive/1.92.zip
	# unzip 1.92.zip
	# cd picard-1.92 && ant -Djava6.home=/usr/lib/jvm/java-6-openjdk-amd64/ -Ddist=/picard -lib lib/ant package-commands
	# picard 1.126
	wget https://github.com/broadinstitute/picard/releases/download/1.126/picard-tools-1.126.zip
	unzip picard-tools-1.126.zip
	sudo cp -r picard-tools-1.126 /picard
