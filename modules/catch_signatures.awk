#! /usr/local/bin/gawk -f

BEGIN {
	if (SZ == "") SZ = 10
	# 默认特征序列大小为 10 bps
}

NR == FNR { # 第一个文件是特征序列总结文件
	# 先复制这个表格到内存中
	sign[FNR]["nm"] = $1
	sign[FNR][1] = $2 ; sign[FNR][2] = $3
	sign[FNR][3] = $4 ; sign[FNR][4] = $5
	sign[FNR]["h1l"] = $7 ; sign[FNR]["h1r"] = $8
	sign[FNR]["h2l"] = $9 ; sign[FNR]["h2r"] = $10
	sign[FNR]["ind"] = $11
	# 这个表格用来存放每个特征序列所在的：
	# nm              - 参考序列名称
	# 1/2/3/4 MH pair - 锚点位置（见下图）
	# h[12][lr]       - 第 1/2 个 MH pair 左/右侧的特征序列
	# ind             - indel 序列（即 MH + internal seq）

	# 再创建一个反向查询表，通过基因组上的锚点反查可能的 MH pair
	# 因为 find_mh 程序在报告多次重复序列中的 MH pair 时，只会输出最左边的 pair
	# 这里通过之前 generate_signature.awk 报告的重复个数，将之后的所有锚点一起加入反查表中
	ids = $4-$2
	for ( i = 0 ; i < $6 ; i++ ) {
		pos[$1, 1, $2+ids*i][FNR] ; pos[$1, 2, $3+ids*i][FNR] ; pos[$1, 3, $4+ids*i][FNR] ; pos[$1, 4, $5+ids*i][FNR]
	}
}

	# position indicator note:
	#           1  2           3  4
	#           |  |           |  |
	#    D site v  v I site    v  v
	# CCTCAGCCAGccgtGTTATAACTTAccgtTTACCAACTACATTTTTTGTAACGAACCAAA
	#           ^ I left clip  |  ^ I right clip
	#              |           ^ D left clip
	#              ^ D right clip

NR != FNR { # 第二个文件为去掉头的 sam 文件，可以用 samtools view 管道传递
	# 备选 reads 为在序列匹配中存在 S（Soft clipping 边缘剪切）I（Insertion 插入）以及 D（Deletion 删除）的 reads
	if ($6 ~ /[SID]/) {
		split($6, seg, "[MSHID]", cig)
		# 首先将 CIGAR 序列分解，seg 中存储的是每个 CIGAR 指示符所代表的子序列长度，cig 存储的是对应的指示符

		# 如果序列存在左剪切
		if ($6 ~ /^[0-9]+S/) {
			# print "read " $1, "Left Clipped" # debug
			if (seg[1] > SZ) {
				# 如果边缘剪切的长度足够生成特征序列，则先在剪切点处先生成这一条 read 的特征序列
				sign_read = substr($10, seg[1]-SZ+1, SZ)
				# 在反查表中查看是否有 MH pair 的锚点在这个剪切点上
				if (length(pos[$3, 1, $4]) > 0) for (i in pos[$3, 1, $4]) {
					# 如果有的话，对于每一个 MH pair 检查其对应的特征序列是否与 read 上的特征序列一致
					if (sign_read == sign[i]["h2l"]) {
						# print "        Found signature #" i, "I LC" # debug
						sign[i]["icount"] ++
					}
				}
				if (length(pos[$3, 3, $4]) > 0) for (i in pos[$3, 3, $4]) {
					if (sign_read == sign[i]["h1l"]) {
						# print "        Found signature #" i, "D LC" # debug
						sign[i]["dcount"] ++
					}
				}
			}
		}
		# 如果存在右剪切
		if ($6 ~ /[0-9]+S$/) {
			# print "read " $1, "Right Clipped" # debug
			if (seg[j-1] > SZ) {
				# 如果边缘剪切的长度足够生成特征序列，则先在剪切点处先生成这一条 read 的特征序列
				sign_read = substr($10, length($10) - seg[j-1] + 1, SZ)
				ep = $4 - 1 # 用于记录 read 在参考序列上的结束位置
				for (j = 1 ; j <= length(cig) ; j++) if (cig[j] == "M" || cig[j] == "D") ep += seg[j]
				# 找到 read 在参考序列上的结束位置，即右剪切的位点
				if (length(pos[$3, 4, ep]) > 0) for (i in pos[$3, 4, ep]) {
					if (seg[j-1] > SZ && sign_read == sign[i]["h1r"]) {
						# print "        Found signature #" i, "I RC" # debug
						sign[i]["icount"] ++
					}
				}
				if (length(pos[$3, 2, ep]) > 0) for (i in pos[$3, 2, ep]) {
					if (seg[j-1] > SZ && sign_read == sign[i]["h2r"]) {
						# print "        Found signature #" i, "D RC" # debug
						sign[i]["dcount"] ++
					}
				}
			}
		}
		# 如果 read 中存在插入序列（必须 indel 位置大小和序列都吻合）
		if ($6 ~ /[0-9]+I/) {
			# print "read " $1, "Ins Inside" # debug
			ins = $4 # 用于记录插入位点在参考序列上的位置
			sp = 1 # 用于记录插入位点在 read 上的位置
			for (j = 1 ; j <= length(cig) ; j++) {
				if (cig[j] == "I") {
					if (length(pos[$3, 1, ins]) > 0) for (i in pos[$3, 1, ins]) { 
						if (sign[i][3]-sign[i][1] == seg[j] && substr($10, sp, seg[j]) == sign[i]["ind"]) {
							# print "        Found signature #" i, "I Inside" # debug
							sign[i]["icount"] ++
						}
					}
				}
				if (cig[j] == "M" || cig[j] == "D") ins += seg[j]
				if (cig[j] == "M" || cig[j] == "I") sp += seg[j]
			}
		}
		# 如果 read 中存在删除序列（只需要 indel 位置大小吻合）
		if ($6 ~ /[0-9]+D/) {
			# print "read " $1, "Del Inside" # debug
			del = $4
			for (j = 1 ; j <= length(cig) ; j++) {
				if (cig[j] == "D") {
					if (length(pos[$3, 1, del]) > 0) for (i in pos[$3, 1, del]) {
						if (sign[i][3]-sign[i][1] == seg[j]) { 
							# print "        Found signature #" i, "D Inside" # debug
							sign[i]["dcount"] ++
						}
					}
				}
				if (cig[j] == "M" || cig[j] == "D") del += seg[j]
			}
		}
	}
}

END {
	OFS = "\t"
	#OFMT = "%8.4f"
	for (i = 1 ; i <= length(sign) ; i++) {
		nm = sign[i]["nm"]
		print nm, sign[i][1], sign[i][2], sign[i][3], sign[i][4], sign[i]["icount"]+0, sign[i]["dcount"]+0
	}
}