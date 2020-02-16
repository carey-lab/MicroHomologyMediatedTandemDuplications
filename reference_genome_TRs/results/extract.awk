BEGIN { print "chr\tpos\tk"}
$1 ~ /^@/  { name = substr($1, 2) }
$1 !~ /^@/ { printf("%s\t%d\t%d\n", name, $1, $3) }