sudo singularity cache clean --force
sudo rm -rf \
    /usr/local/libexec/singularity \
    /usr/local/var/singularity \
    /usr/local/etc/singularity \
    /usr/local/bin/singularity \
    /usr/local/bin/run-singularity \
    /usr/local/etc/bash_completion.d/singularity

sudo rm -r /usr/local/go/src/github.com