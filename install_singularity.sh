# check whether singularity is installed, if not install it
#if sudo apt-get -qq install "singularity-container"; then
  #echo "Singularity already installed"
#else
 # echo "Installing singularity"
  #sudo wget -O- http://neuro.debian.net/lists/focal.de-fzj.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
  #sudo sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
  #sudo apt-get update
  #sudo apt-get install -y singularity-container
#fi

VERSION=3.5.3 # singularity version

# Install dependencies
sudo apt-get update && sudo apt-get install -y \ build-essential \ libssl-dev \ uuid-dev \ libgpgme11-dev \
    squashfs-tools \ libseccomp-dev \ pkg-config \ golang

# Download and compile singularity 3.5.3
sudo mkdir -p /usr/local/go/src/github.com/sylabs && \
  cd /usr/local/go/src/github.com/sylabs && \
  sudo wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
  sudo tar -xzf singularity-${VERSION}.tar.gz && \
  cd ./singularity && \
  sudo ./mconfig && \
  cd builddir && \
  sudo make && \
  sudo make install && \
  cd /usr/local/go/src/github.com/sylabs && \
  sudo rm singularity-${VERSION}.tar.gz