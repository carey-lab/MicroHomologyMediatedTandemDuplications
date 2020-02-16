#! /usr/local/bin/gawk -f
BEGIN {
	OFS = "\t"
	OFMT = "%.2f"
	print "Sequence Name\tH1-left\tH1-right\tH2-left\tH2-right\tRead-Depth\tDup-Read\tin 1M\tCol-read\tin 1M"
}

NR == FNR {
	cov[$1, $2] = $3
}

NR != FNR {
	meancov = int(0.25 * (cov[$1, $2] + cov[$1, $3] + cov[$1, $4] + cov[$1, $5]))
	if (meancov == 0) {
		print $1, $2, $3, $4, $5,     "-", "-",                    "-", "-",                    "-"
	}
	else if ($6 == 0 && $7 == 0) {
		print $1, $2, $3, $4, $5, meancov, "-",                    "-", "-",                    "-"
	}
	else if ($6 == 0) {
		print $1, $2, $3, $4, $5, meancov, "-",                    "-",  $7, $7 / meancov * 1000000
	}
	else if ($7 == 0) {
		print $1, $2, $3, $4, $5, meancov,  $6, $6 / meancov * 1000000, "-",                    "-"
	}
	else {
		print $1, $2, $3, $4, $5, meancov,  $6, $6 / meancov * 1000000,  $7, $7 / meancov * 1000000
	}
}