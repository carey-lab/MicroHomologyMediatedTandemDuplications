initGenomeFromFasta <- function(fa) {
    ## data fields
    fa <- readLines(fa)             ## read the fasta file in all
    ci <- grep("^>", fa)            ## ci = chromasome index
    lw <- nchar(fa[ci + 1])         ## lw = line width in base
    cn <- substring(sapply(strsplit(fa[ci], "\ "), `[[`, 1), 2)
    ln <- lw * (c(ci[-1], length(fa) + 1) - ci - 2) + nchar( fa[c(ci[-1], length(fa) + 1) - 1] )
    content <- list()
    bases <- c("A", "C", "G", "T")
    for (i in 1:length(ci)) {
        s <- unlist(strsplit(fa[(ci[i] + 1):(c(ci[-1], length(fa)+1)[i] - 1)], ""))
        content[[i]] <- c(0, 0, 0, 0) ; names(content[[i]]) <- bases
        for (b in bases) {
            content[[i]][b] <- sum(s == b)
        }
        content[[i]] <- content[[i]] / sum(content[[i]])
    }
    rm(bases, i, s, b)
    # cn = chromasome names
    names(ci) <- cn
    names(lw) <- cn
    names(ln) <- cn
    names(content) <- cn
    # util functions
    seqdiff <- function(s1, s2) {
        # different nucleotide location
        n <- nchar(s1[1])
        if ( nchar(s2[1]) < n ) n <- nchar(s2[1])
        s1 <- substr(s1[1], 1, n)
        s2 <- substr(s2[1], 1, n)
        s <- strsplit(c(s1, s2), "")
        return( which(s[[1]] != s[[2]]) )
    }
    hdist <- function(s1, s2) {
        ## Hamming Distance
        return( length(seqdiff(s1, s2)) )
    }
    extract <- function(chr, start, end = NULL) {
        if (start < 1 || start > ln[chr]) se <- -1
        else if (!is.null(end)) {
            if (end < start || end > ln[chr]) se <- -1
            else {
                # extract sequence from [chr]omasome:[start]-[end]
                se <- paste(fa[(ci[chr] + (start-1) %/% lw[chr] + 1):(ci[chr] + (end-1) %/% lw[chr] + 1)],
                            collapse = "")
                se <- substr(se, (start-1) %% lw[chr] + 1, end - lw[chr] * ((start-1) %/% lw[chr]))
            }
        } else {
            # extract nucleotide from [chr]omosome:[start]-[start]
            se <- substr(fa[(ci[chr] + (start-1) %/% lw[chr] + 1)],
                         (start-1) %% lw[chr] + 1, (start-1) %% lw[chr] + 1)
        }
        if (grepl("N", se)) se <- -1
        return(se)
    }
    ## methods
    list(
        summary = function() {
            for (i in 1:length(cn)){
                cat(cn[i], rep(" ", max(nchar(cn)) - nchar(cn[i])), " (", ln[i] ,") : ", fa[ci[i] + 1],
                    "...\n", sep = "")
            }
        },
        getLength = function(ch) return(ln[ch]),
        getContent = function(ch) return(content[ch]),
        extract = extract
    )
}
