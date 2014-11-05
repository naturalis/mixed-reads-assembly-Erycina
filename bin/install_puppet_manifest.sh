#!/bin/sh
#
# Cloudinit script for deploying mrbayes
#
set -e -x

# Git repository to clone
puppet_source=https://github.com/naturalis/puppet.git

#
# Get latest puppet version
#

# Debian like
if [ -f /usr/bin/dpkg ]
then
  wget http://apt.puppetlabs.com/puppetlabs-release-stable.deb
  dpkg -i puppetlabs-release-stable.deb
  apt-get --yes --quiet update
  apt-get --yes -o Dpkg::Options::="--force-confold" --quiet install git puppet-common ruby1.9.1 libaugeas-ruby
fi

#
# Move original puppet directory
#
if [ -d "/etc/puppet.orig" ]; then
  rm -rf /etc/puppet.orig
fi
mv /etc/puppet /etc/puppet.orig

#
# Fetch puppet configuration from public git repository.
#
env GIT_SSL_NO_VERIFY=true git clone --recursive $puppet_source /etc/puppet

#
# Copy meta data to hiera backend directory
#
# cloud-init json file
if [ -f /meta.js  ]; then
   cp /meta.js /etc/puppet/hieradata/cloud-init.json
fi
# user-data yaml file
if [ -f /user-data.yaml  ]; then
   cp /user-data.yaml /etc/puppet/hieradata/user-data.yaml
fi

puppet apply /etc/puppet/manifests/mrbayes.pp
