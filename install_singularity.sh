# INSTALL SINGULARITY
VERSION=3.5.3 # singularity version
v="$(singularity --version)"

if $v='singularity version 3.5.3'; then
  echo "singularity v3.5.3 already installed"
else
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
fi
