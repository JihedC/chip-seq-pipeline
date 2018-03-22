SHELL=/bin/bash -e

all:
	# Trust the signing key for this repo. Reference: http://cran.rstudio.com/bin/linux/ubuntu/README.html
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
	sudo rm -f /etc/apt/apt.conf.d/99dnanexus
	# Selected Kent utilities v302.1.0
	wget https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedClip
	wget https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedGraphToBigWig
	wget https://github.com/ENCODE-DCC/kentUtils/raw/v302.1.0/bin/linux.x86_64/bedToBigBed
	sudo chmod a+x bedClip bedGraphToBigWig bedToBigBed
	sudo cp bedClip bedGraphToBigWig bedToBigBed /usr/local/bin
