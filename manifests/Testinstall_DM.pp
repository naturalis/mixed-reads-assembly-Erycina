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
	"mk-configure":			ensure => latest;
 	"bzip2":				ensure => latest;
 	"wget":					ensure => latest;
 	"tar":					ensure => latest;
	"zlib1g":				ensure => latest;
	"bwa":					ensure => absent;
	"python":				ensure => latest;
	"python-biopython":		ensure => latest;
	"blast2":				ensure => latest;
	"ncbi-blast+":			ensure => latest;
	"samtools":				ensure => latest;
}

#$username = "dennis"
#$erycina_dir	= "/home/${username}/mixed-reads-assembly-Erycina2"

#file {
#	$erycina_dir:
#		ensure 	=> directory,
#		group 	=> $username,
#		owner 	=> $username,
#		recurse => true;
#}

# command line tasks

exec {
	# install bwa
	"download_bwa":
		command => "wget http://sourceforge.net/projects/bio-bwa/files/latest/download?source=files -O bwa.tar.bz2",
		cwd     => "/home/dennis/Downloads/",
		creates => "/home/dennis/Downloads/bwa.tar.bz2",
		require => Package[ 'wget', 'tar' ];

	"unzip_bwa":
		command => "tar -jxf bwa.tar.bz2",
		cwd     => "/home/dennis/Downloads/",
		creates => "/home/dennis/Downloads/bwa-0.7.10",
		require => Exec['download_bwa'];

	"install_bwa":
		command => "make",
		cwd     => "/home/dennis/Downloads/bwa-0.7.10",
		creates => "/usr/bin/bwa",
		require => Exec['unzip_bwa'];

	"something":
		command => "cp /home/dennis/Downloads/bwa-0.7.10/bwa /usr/bin/",
		creates => "/usr/bin/bwa",
		require => Exec['install_bwa'];
}








