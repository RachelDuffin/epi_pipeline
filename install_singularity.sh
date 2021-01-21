# check whether singularity is installed, if not install it
if sudo apt-get -qq install "singularity-container"; then
  echo "Singularity already installed"
else
  echo "Installing singularity"
  sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
  sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
  sudo apt-get update
  sudo apt-get install -y singularity-container
fi