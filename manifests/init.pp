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
 	"bzip2":					ensure => installed;
 	"wget":						ensure => installed;
 	"tar":						ensure => installed;
	"zlib1g-dev":					ensure => installed;
	"python":					ensure => installed;
	"python-biopython":				ensure => installed;
	"blast2":					ensure => installed;
	"ncbi-blast+":					ensure => installed;
	"samtools":					ensure => installed;


}

# set default paths for data, scripts and source code
$username 	= "janwillem"
$erycina_dir	= "/home/${username}/mixed-reads-assembly-Erycina"  ##"/home/${id}/mixed-reads-assembly-Erycina"
$bin_dir	= "${erycina_dir}/bin"
$data_dir	= "${erycina_dir}/data"
$doc_dir	= "${erycina_dir}/doc"
$doc_paper_dir	= "${doc_dir}/paper"
$manifests_dir	= "${erycina_dir}/manifests"
$results_dir	= "${erycina_dir}/results"
$src_dir	= "${erycina_dir}/src"



# create links for executables and data directories
file {
	$erycina_dir:
		ensure  => directory,
		group   => $username,
		owner   => $username,
		recurse => true;
	$bin_dir:
		ensure  => directory,
		group   => $username,
		owner   => $username,
		recurse => true;
	$data_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
		recurse => true;
	$doc_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
    		recurse => true;
    	$doc_paper_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;
    	$manifests_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
    		recurse => true;
	$results_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;
    	$src_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;
  	"bwa_link":
		path    => "${bin_dir}/bwa",
		ensure  => link,
		target  => "${bin_dir}/bwa-0.7.10",
		require => Exec["unzip_bwa"];
	
}

# command line tasks
exec {
	# install bwa
	"download_bwa":
		command => "wget http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.10.tar.bz2 -O bwa-0.7.10.tar.bz2",
		cwd     => $bin_dir,
		creates => "${bin_dir}/bwa-0.7.10.tar.bz2",
		require => Package[ 'wget', 'tar' ];
	"unzip_bwa":
		command => "tar -jxvf bwa-0.7.10.tar.bz2",
		cwd     => $bin_dir,
		creates => "${bin_dir}/bwa-0.7.10",
		require => Exec['download_bwa'];
	"install_bwa":
		command => "make",
		cwd     => "${bin_dir}/bwa-0.7.10",
		creates => "${bin_dir}/bwa",
		require => Exec['unzip_bwa'];

}
