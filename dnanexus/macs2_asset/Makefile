SHELL=/bin/bash -e

all:
    # Trust the signing key for this repo. Reference: http://cran.rstudio.com/bin/linux/ubuntu/README.html
    sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
    sudo rm -f /etc/apt/apt.conf.d/99dnanexus
    sudo apt-get -y install python-numpy=1:1.11.0-1ubuntu1
    sudo pip install -U --no-deps MACS2==2.1.1.20160309
