# Load necessary packages.
library(tidyverse)


# Load in omega/count tables and factor by gene type.
omega_table_all = read.csv(
    "/Users/davidpeede/Dropbox/Github/santomea_natural_history/data/san_yak_omega_table.csv"
)
omega_table_all$type = factor(
    omega_table_all$type,
    levels = c("Obp", "Csp", "Gr", "Or", "genome"),
    )

mk_table_all = read.csv(
    "/Users/davidpeede/Dropbox/Github/santomea_natural_history/data/san_yak_mk_table.csv"
    )
mk_table_all$type = factor(
    mk_table_all$type,
    levels = c("Obp", "Csp", "Gr", "Or", "genome"),
)


# Define functions.
neutrality_index = function(DS, DA, PS, PA) {
    num = (DS*PA)
    den = (DA*PS)
    NI = (num/den)
    return(NI)
}

ni_sig = function(DS, DA, PS, PA) {
    mat = matrix(c(DS, PS, DA, PA), ncol=2)
    pval = fisher.test(mat, workspace=1e9)$p.value
    return(pval)
}

alpha_bar = function(DS, DA, PS, PA) {
    num = sum((DS*PA)/(PS+DS))
    den = sum((DA*PS)/(PS+DS))
    alpha = 1-(num/den)
    DS_tot = sum(DS)
    DA_tot = sum(DA)
    PS_tot = sum(PS)
    PA_tot = sum(PA)
    mat = matrix(c(DS_tot, PS_tot, DA_tot, PA_tot), ncol=2)
    pval = fisher.test(mat, workspace=1e9)$p.value
    return(alpha, pval)
}

DoS = function(DS, DA, PS, PA) {
    D = (DA)/(DA+DS)
    P = (PA)/(PA+PS)
    DOS = (D-P)
    return(DOS)
}


# Let's look at evoluionary rates first.
# We begin by looking if any specific gene family has elevated omega values
# compared to genomic background.
# To do so we first split the omega data set into san and yak specifc data sets
# removing undefined omega values.
omega_table_san = omega_table_all[!is.na(omega_table_all$omega_san), ]
omega_table_yak = omega_table_all[!is.na(omega_table_all$omega_yak), ]

# First we calculate the number of defined observations per gene family.
genome_total_san = nrow(filter(omega_table_san, type == "genome"))
obp_total_san = nrow(filter(omega_table_san, type == "Obp"))
csp_total_san = nrow(filter(omega_table_san, type == "Csp"))
gr_total_san = nrow(filter(omega_table_san, type == "Gr"))
or_total_san = nrow(filter(omega_table_san, type == "Or"))

genome_total_yak = nrow(filter(omega_table_yak, type == "genome"))
obp_total_yak = nrow(filter(omega_table_yak, type == "Obp"))
csp_total_yak = nrow(filter(omega_table_yak, type == "Csp"))
gr_total_yak = nrow(filter(omega_table_yak, type == "Gr"))
or_total_yak = nrow(filter(omega_table_yak, type == "Or"))

# Next we calculate the number of observations that have an omega value of 1 or
# larger.
genome_pos_san = nrow(filter(omega_table_san, type == "genome" & omega_san >= 1))
obp_pos_san = nrow(filter(omega_table_san, type == "Obp" & omega_san >= 1))
csp_pos_san = nrow(filter(omega_table_san, type == "Csp" & omega_san >= 1))
gr_pos_san = nrow(filter(omega_table_san, type == "Gr" & omega_san >= 1))
or_pos_san = nrow(filter(omega_table_san, type == "Or" & omega_san >= 1))

genome_pos_yak = nrow(filter(omega_table_yak, type == "genome" & omega_yak >= 1))
obp_pos_yak = nrow(filter(omega_table_yak, type == "Obp" & omega_yak >= 1))
csp_pos_yak = nrow(filter(omega_table_yak, type == "Csp" & omega_yak >= 1))
gr_pos_yak = nrow(filter(omega_table_yak, type == "Gr" & omega_yak >= 1))
or_pos_yak = nrow(filter(omega_table_yak, type == "Or" & omega_yak >= 1))

# Next we asses signficance within species for each gene family using a fisher
# exact test.
obp_omega_fet_pval_san = fisher.test(
    matrix(c(obp_pos_san, obp_total_san, genome_pos_san, genome_total_san),
           ncol=2),
)$p.value
csp_omega_fet_pval_san = fisher.test(
    matrix(c(csp_pos_san, csp_total_san, genome_pos_san, genome_total_san),
           ncol=2),
)$p.value
gr_omega_fet_pval_san = fisher.test(
    matrix(c(gr_pos_san, gr_total_san, genome_pos_san, genome_total_san),
           ncol=2),
)$p.value
or_omega_fet_pval_san = fisher.test(
    matrix(c(or_pos_san, or_total_san, genome_pos_san, genome_total_san),
           ncol=2),
)$p.value

obp_omega_fet_pval_yak = fisher.test(
    matrix(c(obp_pos_yak, obp_total_yak, genome_pos_yak, genome_total_yak),
           ncol=2),
)$p.value
csp_omega_fet_pval_yak = fisher.test(
    matrix(c(csp_pos_yak, csp_total_yak, genome_pos_yak, genome_total_yak),
           ncol=2),
)$p.value
gr_omega_fet_pval_yak = fisher.test(
    matrix(c(gr_pos_yak, gr_total_yak, genome_pos_yak, genome_total_yak),
           ncol=2),
)$p.value
or_omega_fet_pval_yak = fisher.test(
    matrix(c(or_pos_yak, or_total_yak, genome_pos_yak, genome_total_yak),
           ncol=2),
)$p.value

# Lastly we asses signficance between species for each gene family using a
# fisher exact test.
obp_omega_fet_pval_between = fisher.test(
    matrix(c(obp_pos_san, obp_total_san, obp_pos_yak, obp_total_yak),
           ncol=2),
)$p.value
csp_omega_fet_pval_between = fisher.test(
    matrix(c(csp_pos_san, csp_total_san, csp_pos_yak, csp_total_yak),
           ncol=2),
)$p.value
gr_omega_fet_pval_between = fisher.test(
    matrix(c(gr_pos_san, gr_total_san, gr_pos_yak, gr_total_yak),
           ncol=2),
)$p.value
or_omega_fet_pval_between = fisher.test(
    matrix(c(or_pos_san, or_total_san, or_pos_yak, or_total_yak),
           ncol=2),
)$p.value
genome_omega_fet_pval_between = fisher.test(
    matrix(c(genome_pos_san, genome_total_san, genome_pos_yak, genome_total_yak),
           ncol=2),
)$p.value

# Next we test for significant differeces between each gene family and the
# genomic background within species.
obp_omega_tt_pval_san = t.test(
    filter(omega_table_san, type == "Obp")$omega_san,
    filter(omega_table_san, type == "genome")$omega_san,
)$p.value
csp_omega_tt_pval_san = t.test(
    filter(omega_table_san, type == "Csp")$omega_san,
    filter(omega_table_san, type == "genome")$omega_san,
)$p.value
gr_omega_tt_pval_san = t.test(
    filter(omega_table_san, type == "Gr")$omega_san,
    filter(omega_table_san, type == "genome")$omega_san,
)$p.value
or_omega_tt_pval_san = t.test(
    filter(omega_table_san, type == "Or")$omega_san,
    filter(omega_table_san, type == "genome")$omega_san,
)$p.value

obp_omega_tt_pval_yak = t.test(
    filter(omega_table_yak, type == "Obp")$omega_yak,
    filter(omega_table_yak, type == "genome")$omega_yak,
)$p.value
csp_omega_tt_pval_yak = t.test(
    filter(omega_table_yak, type == "Csp")$omega_yak,
    filter(omega_table_yak, type == "genome")$omega_yak,
)$p.value
gr_omega_tt_pval_yak = t.test(
    filter(omega_table_yak, type == "Gr")$omega_yak,
    filter(omega_table_yak, type == "genome")$omega_yak,
)$p.value
or_omega_tt_pval_yak = t.test(
    filter(omega_table_yak, type == "Or")$omega_yak,
    filter(omega_table_yak, type == "genome")$omega_yak,
)$p.value

# Lastly we test for significant differeces between each gene family between
# species.
obp_omega_tt_pval_between = t.test(
    filter(omega_table_san, type == "Obp")$omega_san,
    filter(omega_table_yak, type == "Obp")$omega_yak,
)$p.value
csp_omega_tt_pval_between = t.test(
    filter(omega_table_san, type == "Csp")$omega_san,
    filter(omega_table_yak, type == "Csp")$omega_yak,
)$p.value
gr_omega_tt_pval_between = t.test(
    filter(omega_table_san, type == "Gr")$omega_san,
    filter(omega_table_yak, type == "Gr")$omega_yak,
)$p.value
or_omega_tt_pval_between = t.test(
    filter(omega_table_san, type == "Or")$omega_san,
    filter(omega_table_yak, type == "Or")$omega_yak,
)$p.value
genome_omega_tt_pval_between = t.test(
    filter(omega_table_san, type == "genome")$omega_san,
    filter(omega_table_yak, type == "genome")$omega_yak,
)$p.value


# Lets look at the MK count data.
# Now lets calculate the neutrality index. First we have to remove any gene
# that has a value of 5 or less in any cell to avoid biased results.
ni_table_san = filter(
    mk_table_all,
    (DA_san >= 5) & (DS_san >= 5) & (PA_san >= 5) & (PS_san >= 5),
)
ni_table_yak = filter(
    mk_table_all,
    (DA_yak >= 5) & (DS_yak >= 5) & (PA_yak >= 5) & (PS_yak >= 5),
)

# Next we calculate the neutrality index per gene and asses significance using
# a Fisher's Exact Test.
ni_table_san$NI = neutrality_index(
    ni_table_san$DS_san, ni_table_san$DA_san,
    ni_table_san$PS_san, ni_table_san$PA_san
)
ni_table_san$pval = ni_sig(
    ni_table_san$DS_san, ni_table_san$DA_san,
    ni_table_san$PS_san, ni_table_san$PA_san
)


















