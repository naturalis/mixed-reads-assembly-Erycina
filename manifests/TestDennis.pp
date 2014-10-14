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
package {
 	"bzip2":					ensure => installed;
 	"wget":						ensure => installed;
 	"tar":						ensure => installed;
	"zlib1g":					ensure => installed;
	"bwa":						ensure => installed;
	"python":					ensure => installed;
	"python-biopython":				ensure => installed;
	"blast2":					ensure => installed;
	"ncbi-blast+":					ensure => installed;
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
}
