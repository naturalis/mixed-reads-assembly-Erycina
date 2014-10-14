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
	"python-biopython":			ensure => installed;
	"blast2":					ensure => installed;
	"ncbi-blast+":				ensure => installed;
}
