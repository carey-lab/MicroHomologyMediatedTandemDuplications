#### Utility for TR pattern and flanking MH identification
#### This utility requires the genome utility clousure

# sequence identity below this threshold will be denied
consensus_thresh <- 0.75
consensus_lower <- 0.5
# consensus_thresh <- 0.3
# consensus_lower <- 0.5
# 4 bases
bases <- c("A", "C", "G", "T")

# this function returns the base freq matrix given a vector of
# sequences with the same length (nchar)
makeFrequencyMatrix <- function(seqs) {
    n <- nchar(seqs[1])
    seqs <- matrix(unlist(strsplit(seqs, "")), n)
    m <- matrix(0, 4, n)
    for (j in 1:n) {
        b <- seqs[j, ]
        for (i in 1:4) {
            m[i, j] <- sum(b == bases[i])
        }
    }
    rownames(m) <- bases
    m / ncol(seqs)
}
# this function returns a consensus sequence given a freq matrix
makeConsensus <- function(fm) {
    use <- paste(apply(fm, 2, function(cl) bases[which(cl == max(cl))[1]]), collapse = "")
    out <- paste(apply(fm, 2, function(cl) {
        if (all(cl < 0.4)) b <- "*" else {
            b <- bases[which(cl > max(cl) - 0.1)]
            if (length(b) > 1) b <- paste0("[", paste(b, collapse = "/"), "]")
            else b <- bases[which(cl == max(cl))[1]]
        }
        b
    } ), collapse = "")
    attr(use, "output") <- out
    use
}
# this function returns the entropy of a given freq matrix
seqEntropy <- function(fm) {
    baserate <- rowSums(fm) / ncol(fm)
    sum(ifelse(baserate == 0, 0, -baserate * log2(baserate)))
}
# this function returns the information content of a given freq matrix
# infoContent <- function(fm) {}
# this function calculate the match between a sequence and a freq matrix
match <- function(seq, fm) {
    seq <- unlist(strsplit(seq, ""))
    s <- 0
    for (i in 1:length(seq)) {
        s <- s + fm[seq[i], (i-1)%%ncol(fm)+1]
    }
    s
}
# this function find if a given subsequence have inner repetitive pattern
findInnerRep <- function(seq) {
    # find the inner repeating pattern of a given seq
    n <- nchar(seq)
    inner <- 1
    for (i in (2:n)[n %% (2:n) == 0]) {
        ss <- substring(seq, (1:i - 1) * n/i + 1, 1:i * n/i)
        fm <- makeFrequencyMatrix(ss)
        if (all(sapply(ss, match, fm = fm) > consensus_thresh * n / i))
            inner <- i
    }
    inner
}

# this is the main function that identifies the TR given the unit info
annotateTR <- function(gu, chr, posOriginal, k) {
    pos <- posOriginal
    seq <- gu$extract(chr, pos, pos+k-1)
    endpos <- pos + k
    
    ir <- 1
    ir <- tryCatch( { findInnerRep(seq) })
    if (ir > 1) seq <- substring(seq, (1:ir - 1) * nchar(seq)/ir + 1, 1:ir * nchar(seq)/ir)
    k <- nchar(seq[1])

    lor <- 0
    score <- numeric()
    repeat {
        e <- gu$extract(chr, pos - (lor+1) * k, pos - lor * k - 1)
        if (e == -1) break
        else {
            s <- match(e, makeFrequencyMatrix(seq))
            if (s < consensus_lower*k) break
            else { score[length(score)+1] <- s ; seq <- c(e, seq) ; lor <- lor+1 }
        }
    }
    lor <- which(score > consensus_thresh*k)
    if (length(lor)) lor <- max(lor) else lor <- 0
    seq <- seq[(length(seq)-ir-lor+1):length(seq)]
    
    ror <- 0
    score <- numeric()
    repeat {
        e <- gu$extract(chr, endpos + ror * k, endpos + (ror+1) * k - 1)
        if (e == -1) break
        else {
            s <- match(e, makeFrequencyMatrix(seq))
            if (s < consensus_lower*k) break
            else { score[length(score)+1] <- s ; seq <- c(seq, e) ; ror <- ror+1 }
        }
    }
    ror <- which(score > consensus_thresh*k)
    if (length(ror)) ror <- max(ror) else ror <- 0
    seq <- seq[1:(lor+ir+ror)]
    
    pos <- pos - lor * k
    rep <- length(seq)
    
    if (rep > 1) {
        repeat {
            e <- gu$extract(chr, pos - 1)
            if (e == -1) break
            else {
                shift <- paste0( c(e, substring(seq[1:(rep-1)], k, k)), substring(seq, 1, k-1) )
                if (match(shift, makeFrequencyMatrix(shift)) >= match(seq, makeFrequencyMatrix(seq))) {
                    pos <- pos-1 ; seq <- shift 
                    e <- gu$extract(chr, pos + k*rep, pos + k*(rep+1) - 1)
                    if (e != -1 && match(e, makeFrequencyMatrix(seq)) >= consensus_thresh*k) {
                        seq <- c(seq, e) ; rep <- rep+1
                    }
                } else break
            }
        }
    }
    
    fl <- 0
    mh <- ""
    counter_to_avoid_loop = 0
    repeat {
        counter_to_avoid_loop = counter_to_avoid_loop + 1 
        if(counter_to_avoid_loop>1000){
            mh = 'XXXXXXXXXXXXXXXXXXXXXX'
            break
        }
        if (fl >= consensus_thresh*k && k-fl < 4) {
            e <- gu$extract(chr, pos + rep * k, pos + (rep+1) * k-1)
            if (e != -1) {
                seq[rep+1] <- gu$extract(chr, pos + rep * k, pos + (rep+1) * k-1)
                rep <- rep + 1 ; fl <- 0 ; mh <- ""
            }
        } else {
            e <- gu$extract(chr, pos + rep * k, pos + rep * k + fl)
            if (e == -1) break
            else {
                if (all(e == substring(seq, 1, fl+1))) {
                    fl <- fl+1 ; mh <- e
                } else break
            }
        }
    }
    if (nchar(mh) < 2) mh <- ""
    
    mat <- makeFrequencyMatrix(seq)
    cons <- makeConsensus(mat)
    entropy <- seqEntropy(mat)
    
    list(chr = chr, posOriginal = posOriginal , pos = pos, rep = rep, k = k, seq = seq,
         cons = attr(cons, "output"), entropy = entropy, mh = mh)
}
# annotateIndelTR <- function(gu, chr, pos, seq, is_insert = F) {
#     endpos <- ifelse(is_insert, pos, pos + nchar(seq))
#     ir <- findInnerRep(seq) # repeat unit count in the indel seq
#     if (ir > 1) seq <- substring(seq, (1:ir - 1) * nchar(seq)/ir + 1, 1:ir * nchar(seq)/ir)
#     k <- nchar(seq[1])
#     lor <- 0 # repeat unit count on the left side of indel site
#     while ({
#         e <- gu$extract(chr, pos - (lor+1) * k + 1, pos - lor * k)
#         if (e != -1) {
#             seq <- c(e, seq)
#             meanEntropy(makeFrequencyMatrix(seq)) < consensus_thresh
#         } else {
#             if (1 < pos - lor * k) {
#                 seq <- c(gu$extract(chr, 1, pos - lor * k), seq)
#             } else seq <- c("", seq)
#             FALSE
#         }
#     }) lor <- lor + 1
#     ror <- 0 # repeat unit count on the right side of indel site
#     while ({
#         e <- gu$extract(chr, endpos + ror * k, endpos + (ror+1) * k - 1)
#         if (e != -1) {
#             seq <- c(seq, e)
#             meanEntropy(makeFrequencyMatrix(seq[-1])) < consensus_thresh
#         } else {
#             if (endpos + ror * k < gu$getLength(chr)) {
#                 seq <- c(seq, gu$extract(chr, endpos + ror * k, gu$getLength(chr)))
#             } else seq <- c(seq, "")
#             FALSE
#         }
#     }) ror <- ror + 1
#     mat <- makeFrequencyMatrix(seq[2:(length(seq)-1)]) # frequency matrix
#     con <- makeConsensus(mat) # consensus sequence
#     lf <- 0 # left flanking homology
#     # while ({
#     #     pb <- mat[substring(seq[1], k-lf, k-lf), k-lf]
#     #     # pb > 0.4 && pb > max(mat[, k-lf]) - 0.1
#     #     pb == 1
#     # }) lf <- lf + 1
#     rf <- 0 # right flanking homology
#     while ( (rf < nchar(seq[length(seq)])) && 
#             (mat[substring(seq[length(seq)], rf+1, rf+1), rf+1] == 1) ) rf <- rf + 1
#     if (lf + rf >= 3 && lf + rf <= min(15, 0.75 * k) ) { # recognize as microhomology
#         if (lf > 0 && rf == 0) mh <- mat[ ,(k-lf+1):k]
#         else if (lf ==0 && rf > 0) mh <- mat[, 1:rf]
#         else mh <- mat[, c((k - lf + 1):k, 1:rf)]
#         mh <- makeConsensus(mh)
#         if (lf != 0) seq[1:(length(seq) - 1)] <- paste(
#             substring(seq[1:(length(seq) - 1)], 1, k-lf), 
#             tolower(substring(seq[1:(length(seq) - 1)], k-lf+1, k)), sep = ""
#         )
#         if (rf != 0) seq[2:length(seq)] <- paste(
#             tolower(substring(seq[2:length(seq)], 1, rf)),
#             substring(seq[2:length(seq)], rf+1, k), sep = ""
#         )
#     } else mh <- ""
#     rnc <- if (is_insert) { c(lor + ror, lor + ir + ror)
#     } else { c(lor + ir + ror, lor + ror) } # repeat number
#     # reformat sequence
#     if (rnc[1] > 0) {
#         seq[(lor+2)] <- paste0(ifelse(is_insert, ">", "<"), seq[(lor+2)])
#         seq[(lor+1+ir)] <- paste0(seq[(lor+1+ir)], ifelse(is_insert, "<", ">"))
#         seq <- paste0(chr, ":",pos - (lor+1) * k, " - ", paste(seq, collapse = " "))
#         repstarts <- pos + c((-lor):(ir + ror - is_insert - 1)) * k
#         repends <- repstarts + k
#         mhstarts <- c(repstarts[1], repends)
#         mhends <- mhstarts + nchar(mh)
#         return(list(repeat_number = rnc, consensus = con, seq = seq, microhomology = mh,
#                     repeat_starts = repstarts, repeat_ends = repends,
#                     mh_starts = mhstarts, mh_ends = mhends))
#     } else return(NULL)
# }
