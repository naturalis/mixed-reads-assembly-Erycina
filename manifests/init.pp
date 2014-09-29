    #JDOUBLE-U
      file { 'ntp.conf':
        path    => '/etc/ntp.conf',
        ensure  => file,
        require => Package['ntp'],
        source  => "puppet:///modules/ntp/${conf_file}",
      }
    }
