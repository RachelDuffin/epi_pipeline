# INSTALL SINGULARITY
VERSION=3.5.3 # singularity version
v="$(singularity --version)"

if $v='singularity version 3.5.3'; then
  echo "singularity v3.5.3 already installed"
else
  # Install checkinstall
  sudo apt-get update && sudo apt-get install checkinstall
  # Install dependencies
  sudo apt-get update && sudo apt-get install -y build-essential libssl-dev uuid-dev libgpgme11-dev squashfs-tools libseccomp-dev pkg-config golang
  # Download and compile singularity 3.5.3
  mkdir -p /usr/local/go/src/github.com/sylabs && \
  cd /usr/local/go/src/github.com/sylabs && \
  wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
  tar -xzf singularity-${VERSION}.tar.gz && \
  cd singularity && \
  ./mconfig && \
  cd builddir && \
  make && \
  sudo checkinstall && \
  cd /usr/local/go/src/github.com/sylabs && \
  rm singularity-${VERSION}.tar.gz
fi
