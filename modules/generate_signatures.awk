#!/usr/local/bin/gawk -f
# usage
# <SZ=[feature size]> reference.fa mh.1.out <mh.2.out ...>

BEGIN {
	if (SZ == "") SZ = 10 # 默认特征序列的大小为 10bps
	printf("Feature size:          %7d\n", SZ) | "cat >& 2"
	print "-------------------------------------------" | "cat >& 2"
	nname = 0;
	iname = 0;
	tot = 0 ; outrange = 0 ; fail = 0
}

NR == FNR { # 首先将 fasta 文件读入内存
	if ($0 ~ "^>") {
		nname++
		sname[nname] = substr($1, 2)
		block = 0
	}
	else {
		if (block == 0) g[sname[nname], "b"] = length($0)
		g[sname[nname], block] = $0
		block += length($0)
		g[sname[nname], "l"] = block
	}
	# 表格 g 中，每个序列名称下：
	# g[name, "b"]       - 此序列的行宽
	# g[name, "l"]       - 此序列的长度
	# g[name, <int> * b] - 第 <int> 行的序列，如 g[name, 120] 如果行宽为 60 那么这个项中存储的是第三行的序列
}

# 这个函数返回在 g 中储存的，序列名称为 name 的序列，从 s 到 e 碱基之间的序列
# 如果超出范围，则返回 -1
function g_extract(name, s, e,
	ss, ee, i, b, r) {
	if (g[name, "l"] == "") return - 1
	if (s < 1) return -1
	if (e > g[name, "l"]) return -1
	if (s > e) return -1
	b = g[name, "b"]
	s--
	e--
	ss = s - (s % b)
	ee = e - (e % b)
	# print name, b, s, ss, e, ee
	if (ss == ee) {
		r = substr(g[name, ss], s % b + 1, e - s + 1)
	}
	else {
		r = substr(g[name, ss], s % b + 1)
		for (i = ss + b ; i < ee ; i += b) {
			r = r g[name, i]
		}
		r = r substr(g[name, ee], 1, e % b + 1)
	}
	return r
}

# 这个函数用来从 g 中储存的，序列名称为 n 的序列中，提取 H1 和 H2 左位点为 s 和 e，MH 长度为 k 的 MH pair 的特征序列
# 这个函数将结果返回到名为 r 的 array 中
function gen_signature(n, s, e, k, r,
	h1l, h1r, h2l, h2r, ind, dis, i, d) {
	# 首先先确定 MH 位置是否有 Tandem Repeats
	d = e - s # 这个是 MH pair 的 indel size
	for (i = 1 ; 1 ; i++) {
		if (g_extract(n, s, e+k-1) != g_extract(n, s+d*i, e+d*i+k-1)) break
	}
	# 然后是特征序列 H1 的 left/right H2 的 left/right
	h1l = g_extract(n, s-SZ, s-1)
	h1r = g_extract(n, s+k, s+k+SZ-1)
	h2l = g_extract(n, e-SZ, e-1)
	h2r = g_extract(n, s+k+i*d, s+k+i*d+SZ-1) # 这里要求的是最后一个 TR 右侧的特征序列
	ind = g_extract(n, s, e-1)
	# 如果有超出序列范围的情况 报错退出
	if (h1l == -1 || h1r == -1 || h2l == -1 || h2r == -1) return -1
	
	r[0] = i
	r[1] = h1l ; r[2] = h1r ; r[3] = h2l ; r[4] = h2r ; r[5] = ind
	return 0
}

NR != FNR {
	if (FNR == 1) iname ++
	# print SZ, SN, $1, $2, $3
	tot ++
	test = gen_signature(sname[iname], $1, $2, $3, sign)
	if (test == 0) {
		print sname[iname], $1, $1+$3-1, $2, $2+$3-1, sign[0], sign[1], sign[2], sign[3], sign[4], sign[5]
		# the output records are:
		# sequence name | I/D clip sites 1 2 3 4 |
		# WT signature 1 | WT signature 2 | I signature | D signature
		# position indicator note:
		#           1  2           3  4
		#           |  |           |  |
		#    D site v  v I site    v  v
		# CCTCAGCCAGccgtGTTATAACTTAccgtTTACCAACTACATTTTTTGTAACGAACCAAA
		#           ^ I left clip  |  ^ I right clip
		#              |           ^ D left clip
		#              ^ D right clip
	}
	else (test == -1) outrange ++
	# else if (test == -2) fail ++
}

END {
	printf("Total MH pairs:        %7d\n", tot) | "cat >& 2"
	printf("Successed:             %7d (%.1f%)\n", tot - outrange - fail, (tot - outrange - fail) / tot * 100)| "cat >& 2"
	printf("Signature out of range:%7d (%.1f%)\n", outrange, outrange / tot * 100) | "cat >& 2"
	# printf("Specificity failed:    %7d (%.1f%)\n", fail, fail / tot * 100) | "cat >& 2"
}
