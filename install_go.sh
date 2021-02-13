# Download Go
GO_VERSION="go1.15.8.linux-amd64.tar.gz"
cd apps
if [ -f "$GO_VERSION" ] ; then
  echo "tar file already downloaded"
else
  wget https://golang.org/dl/${GO_VERSION}
fi

# Install Go
sudo tar -C /usr/local -zxvf ${GO_VERSION}

# Set up environment for go
export GOPATH=$HOME/work
export PATH=$PATH:/usr/local/go/bin:$GOPATH/bin


source ~/.bash_profile


go version
# now reboot your computer
