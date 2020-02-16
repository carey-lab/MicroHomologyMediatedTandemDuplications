## prepare genome utils
source("genome_util_lite.R")
gu_sp <- initGenomeFromFasta("data/s_pombe.fa")
#gu_sc <- initGenomeFromFasta("data/s_cerevisiae.fa")

## read TR utils
source("code/tr_util.R")

## prepare data for pombe
# read TRF output summarized table
pombe_orig <- read.delim("results/s_pombe.out")
# read base-wise genomic feature mapping
pombe_feature <- list(I   = readLines("~/lab-projects/li-yuze/MHp-finding/results/features.new.I.out"),
                      II  = readLines("~/lab-projects/li-yuze/MHp-finding/results/features.new.II.out"),
                      III = readLines("~/lab-projects/li-yuze/MHp-finding/results/features.new.III.out"))
# revise the result using homemade TR util
pombe <- pombe_orig
for (i in 1:nrow(pombe)) {
    summary <- annotateTR(gu_sp, pombe[i, "chr"], pombe[i, "pos"], pombe[i, "k"])
    pombe$pos[i] <- summary$pos; pombe$k[i] <- summary$k; pombe$rep[i] <- summary$rep;
    pombe$consensus[i] <- summary$cons ; pombe$entropy[i] <- summary$entropy
    pombe$mh[i] <- summary$mh != "" ; pombe$mh_seq[i] <- summary$mh
    f <- pombe_feature[[summary$chr]][summary$pos:(summary$pos + summary$k * summary$rep - 1)]
    pombe$region[i] <- names(table(f))[1]
}
pombe <- unique(pombe)
pombe_filter <- data.frame(k = pombe$k < 10, e = pombe$entropy < 1.5, r = pombe$rep < 2)
pombe_use <- pombe[!(pombe_filter$k | pombe_filter$e | pombe_filter$r), ]

regions <- c("intergenic", "intron", "CDS", "5-UTR", "3-UTR", "ncRNA", "pseudogene")
off <- 0.025*nrow(pombe_use)

plot_tab <- function(table, labx, laby, vlab = F, title = "", titline = 1.2, layoutwidth = 1) {
    off <- off / layoutwidth
    pl <- list(
        x1 = c(0, cumsum(rowSums(table[1:(nrow(table)-1), ,drop=FALSE]))) + (1:nrow(table)-1) * off,
        x2 = cumsum(rowSums(table)) + (1:nrow(table)-1) * off,
        y = table[, 2] / rowSums(table)
    )
    plot(c(0, pl$x2[length(pl$x2)]), c(0, 1.025),
         ann = F, axes = F, type = "n")
    for (i in 1:length(pl$x1)) {
        polygon(c(pl$x1[i], pl$x2[i], pl$x2[i], pl$x1[i]),
                c(0, 0, pl$y[i], pl$y[i]), col = "grey80")
        polygon(c(pl$x1[i], pl$x2[i], pl$x2[i], pl$x1[i]),
                c(pl$y[i], pl$y[i], 1, 1) + 0.025)
        if (pl$x2[i]-pl$x1[i] > 4) {
            text(mean(c(pl$x1[i], pl$x2[i])), pl$y[i]-0.04 , table[i, 2],
                 cex = 0.75, xpd = T)
            text(mean(c(pl$x1[i], pl$x2[i])), pl$y[i]+0.065, table[i, 1],
                 cex = 0.75, xpd = T)
        } else {
            text(pl$x2[i]+0.5*off, pl$y[i]-0.04 , table[i, 2],
                 cex = 0.75, xpd = T)
            text(pl$x2[i]+0.5*off, pl$y[i]+0.065, table[i, 1],
                 cex = 0.75, xpd = T)
        }
    }
    if (vlab) axis(1, at = (pl$x1+pl$x2)/2, tick = F, labels = labx, line = -0.6, las = 2)
    else      axis(1, at = (pl$x1+pl$x2)/2, tick = F, labels = labx, line = -0.6)
    axis(2, at = c(0.1, 1.025-0.1), tick = F, labels = laby, line = -1)
    title(xlab = title, line = titline)
}


opar <- par(no.readonly = T)

png("plot/figure_main.png", width = 7.5, height = 2.7, unit = "in", res = 300)
par(mar = c(6.1, 2.1, 0.6, 0.3), oma = c(0, 0, 1, 0))
layout(matrix(c(1,2,3), nrow = 1), widths = c(5, 5, 3))

# summary tables for plotting panel A
pombe_rep <- cut(pombe_use$rep, breaks = c(2, 3, 4, 7, Inf), right = F, labels = c("2", "3", "4-6", ">6"))
pombe_mh <- pombe_use$mh
tab_rep_mh <- table(pombe_rep, pombe_mh)
plot_tab(tab_rep_mh, labx = levels(pombe_rep), laby = c("with MH", "no MH"), title = "Repeat Number")


# summary tables for plotting panel B
pombe_region <- pombe_use$region
tab_reg_mh <- table(pombe_region, pombe_mh)
plot_tab(tab_reg_mh, labx = regions, laby = c("with MH", "no MH"), vlab = T)


# summary tables for panel C
pombe_m3 <- pombe_use$k %% 3 == 0
tab_mh_m3 <- table(pombe_mh, pombe_m3)
plot_tab(tab_mh_m3, labx = c("no MH", "with MH"), laby = c("unit size 3-divisible", "not 3-divisible"), layoutwidth = 3/(5+3))

dev.off()



png("plot/figure_sup.png", width = 7.5, height = 2.7, unit = "in", res = 300)
par(mar = c(6.1, 2.1, 0.6, 0.3), oma = c(0, 0, 1, 0))
layout(matrix(c(1,2,3), nrow = 1), widths = c(5, 5, 3))

# summary table for sup
tab_reg_m3 <- table(pombe_region, pombe_m3)
plot_tab(tab_reg_m3, labx = regions, laby = c("unit size 3-divisible", "not 3-divisible"), vlab = T)

dev.off()

