SHELL=/bin/bash -e

all:
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	sudo apt-get -y install openjdk-6-jdk=6b41-1.13.13-0ubuntu0.14.04.1
	wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
	unzip Trimmomatic-0.36.zip
	sudo cp -r Trimmomatic-0.36/ /trimmomatic
