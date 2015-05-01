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
 	"bzip2":					ensure => latest;
 	"wget":						ensure => latest;
 	"tar":						ensure => latest;
	"zlib1g-dev":					ensure => latest;
	"python":					ensure => latest;
	"python-biopython":				ensure => latest;
	"blast2":					ensure => latest;
	"ncbi-blast+":					ensure => latest;
	"samtools":					ensure => latest;
	"make":						ensure => latest;
	"seqtk":					ensure => latest;

}


# set default paths for data, scripts and source code
$username 			= "ubuntu"
$erycina_dir			= "/home/${username}/mixed-reads-assembly-Erycina"  ##"/home/${id}/mixed-reads-assembly-Erycina"
$bin_dir			= "${erycina_dir}/bin"
$data_dir			= "${erycina_dir}/data"
$doc_dir			= "${erycina_dir}/doc"
$doc_paper_dir			= "${doc_dir}/paper"
$results_dir			= "${erycina_dir}/results"
$results_assembly_dir		= "${results_dir}/assembly"
$results_local_blast_dir	= "${result_dir}/local_blast"
$src_dir			= "${erycina_dir}/src"


# create links for executables directories
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
	$results_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;
	$results_assembly_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
    		recurse => true;
	$results_local_blast_dir:
		ensure  => directory,
    		group   => $username,
    		owner   => $username,
    		recurse => true;
    	$src_dir:
		ensure  => directory,
	  	group   => $username,
    		owner   => $username,
    		recurse => true;
#  	"bwa_link":
#		path    => "${bin_dir}/bwa",
#		ensure  => link,
#		target  => "${bin_dir}/bwa-0.7.10",
#		require => Exec["unzip_bwa"];
	
}

# command line tasks
exec {
	# install bwa
	"download_bwa":
		command => "wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download -O bwa.tar.bz2",
		cwd     => $bin_dir,
		creates => "${bin_dir}/bwa.tar.bz2",
		require => Package[ 'wget', 'tar' ];
	"unzip_bwa":
		command => "tar -jxvf bwa.tar.bz2",
		cwd     => $bin_dir,
		creates => "${bin_dir}/bwa-0.7.10",
		require => Exec['download_bwa'];
	"make_bwa":
		command => "make",
		cwd     => "${bin_dir}/bwa-0.7.10",
		creates => "${bin_dir}/bwa",
		require => Exec['unzip_bwa'];
	"move_exe":
		command => "cp ${bin_dir}/bwa-0.7.10/bwa /usr/bin/",
		creates => "/usr/bin/bwa",
		require => Exec['make_bwa'];

}
