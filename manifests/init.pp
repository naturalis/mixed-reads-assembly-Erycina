# update the $PATH environment variable for the Exec tasks.
 Exec {
 	path => [
 		"/usr/local/sbin",
 				"/usr/local/bin",
 		"/usr/sbin",
 		"/usr/bin",
 		"/sbin",
 		"/bin",
 	]
 }

 # Installing required packages.
 package {
 			"bunzip2":                  ensure => installed;
 	"wget":                     ensure => installed;
 	"tar":                      ensure => installed;
    	"bwa":                      ensure => installed;
    	"zlib-devel":		    ensure => installed;

}

# set default paths for storing data, tools and source code
$username = "hettling"
$erycina_dir	= "/home/${username}/assembly_Erycina"  ##"/home/${id}/assembly_Erycina"
$tools_dir	= "${erycina_dir}/tools"
$tools_bin_dir	= "${tools_dir}/bin"
$src_dir	= "${erycina_dir}/src"
$data_dir	= "${erycina_dir}/data"

# create links for executables and data directories
file {
	$erycina_dir:
		ensure  => directory,
		group   => $username,
		owner   => $username,
		recurse => true;
	$data_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
		recurse => true;
	$src_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
    		recurse => true;
	$tools_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;

  	$tools_bin_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;

  	"muscle_link":
		path    => "${tools_bin_dir}/muscle",
		ensure  => link,
		target  => "${tools_bin_dir}/muscle3.8.31_i86linux64",
		require => Exec["unzip_muscle"];
	"consense_link":
		path    => "/usr/local/bin/consense",
		ensure  => link,
		target  => "${tools_dir}/phylip-3.696/src/consense",
		require => Exec["make_phylip"];
	"examl_link":
		path    => "/usr/local/bin/examl",
		ensure  => link,
		target  => "${tools_dir}/ExaML/examl/examl",
		require => Exec["compile_examl"];
    	"exabayes_link":
		path    => "/usr/local/bin/exabayes",
		ensure  => link,
		target  => "${tools_dir}/exabayes-1.2.1/bin/exabayes",
		require => Exec["compile_exabayes"];
        "treepl_link":
	    	path    => "/usr/local/bin/treePL",
		ensure  => link,
		target  => "${tools_dir}/treePL/src/treePL",
		require => Exec["compile_treepl"];
	"inparanoid_dir":
		path    => "${data_dir}/inparanoid",
		ensure  => directory;
}

# command line tasks
exec {
  # add bin directory for all required tools to PATH
  "make_bindir_sh":
    command => "echo 'export PATH=\$PATH:${tools_bin_dir}' > supersmart-tools-bin.sh",
    cwd     => "/etc/profile.d",
    creates => "/etc/profile.d/supersmart-tools-bin.sh",
    require => Exec[ 'clone_bio_phylo' ];
  "make_bindir_csh":
    command => "echo 'setenv PATH \$PATH:${tools_bin_dir}' > supersmart-tools-bin.csh",
    cwd     => "/etc/profile.d",
    creates => "/etc/profile.d/supersmart-tools-bin.csh",
    require => File[ $tools_bin_dir ];
	# install bwa
	"download_bwa":
		command => "wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.10.tar.bz2 -O bwa-0.7.10.tar.bz2",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bwa-0.7.10.tar.bz2",
		require => Package[ 'wget', 'tar' ];
        "bunzip_bwa":
		command => "bunzip2 bwa-0.7.10.tar.bz2",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bwa-0.7.10.tar.bz2",
		require => Package[ 'bunzip2' ];
	"unzip_bwa":
		command => "tar -xvf bwa-0.7.10.tar",
		cwd     => $tools_dir,
		creates => "${tools_dir}/bwa-0.7.10.tar",
		require => Exec['download_bwa'];
	"install_bwa":
		command => "${tools_dir}/bwa-0.7.10/configure --prefix=/usr --disable-dlopen && make install",
		creates => "/usr/bin/bwa",
		cwd     => "${tools_dir}/bwa-0.7.10",
		timeout => 0,
		require => Exec['unzip_bwa'];

}
