#' @title Annotate mutations with their sequence-motif context
#'
#' @description Labels each single-nucleotide variant by whether it falls at
#' the target base of a sequence motif (e.g. the C of the AID hotspot WRCY),
#' elsewhere within the motif, or outside it. Both strands are checked.
#'
#' @details For each SNV the reference sequence around the mutated base is
#' retrieved (once per distinct position). A variant is labelled:
#' \describe{
#'   \item{SITE}{the mutated base is the `index` base of the motif on the
#'   forward strand, or the corresponding base of its reverse complement
#'   (e.g. the G of RGYW for WRCY);}
#'   \item{MOTIF}{the mutated base lies inside an occurrence of the motif
#'   (either strand) but not at the indexed base;}
#'   \item{FALSE}{neither.}
#' }
#' Indels and multi-base substitutions are not assessed and get `NA`.
#'
#' With `aid_alleles = TRUE`, a SITE also requires the substitution expected
#' from AID/error-prone repair, read on the motif's strand: C>T or C>G at a C
#' target, A>G or A>T at an A target (e.g. the WA motif of polymerase eta).
#' Substitutions at the site with other alleles are then labelled MOTIF.
#'
#' The reference comes from `fastaPath` (an indexed FASTA) when given,
#' otherwise from the BSgenome package for `genome_build`
#' (BSgenome.Hsapiens.UCSC.hg19 for grch37, BSgenome.Hsapiens.UCSC.hg38 for
#' hg38). Chromosome names with or without a "chr" prefix are matched to the
#' reference's naming automatically.
#'
#' @param maf MAF data frame (required columns: Chromosome, Start_Position,
#'   Reference_Allele, Tumor_Seq_Allele2).
#' @param motif Motif in IUPAC code (default "WRCY"). A vector annotates
#'   several motifs at once, adding one column per motif.
#' @param index Position of the target base within the motif (default 3, the
#'   C of WRCY). Recycled to the length of `motif`.
#' @param genome_build "grch37" or "hg38". Taken from `maf` when it is a
#'   maf_data object.
#' @param fastaPath Optional path to an indexed FASTA of the matching build.
#' @param aid_alleles If TRUE, require an AID-type substitution for SITE (see
#'   Details). Default FALSE.
#' @param motif_strand If TRUE, add a `<motif>_strand` column giving the strand
#'   on which the motif was found for SITE calls: "+" (the reference base is
#'   the motif's target) or "-" (it is the target's complement on the reverse
#'   strand). Default FALSE.
#'
#' @return `maf` with a `seq` column (reference sequence centred on the
#'   variant, spanning the first motif's length either side) and one column per
#'   motif, named after the motif, with values "SITE", "MOTIF", "FALSE" or NA.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' my_maf <- GAMBLR.open::get_coding_ssm() %>%
#'   dplyr::filter(Hugo_Symbol=="BCL2")
#'
#' annotated = annotate_ssm_context(maf = my_maf,
#'                                  motif = c("WRCY", "AGCT", "WA"),
#'                                  index = c(3, 3, 2),
#'                                  genome_build = "grch37",
#'                                  aid_alleles = TRUE)
#' table(annotated$WRCY, useNA = "ifany")
#'}
annotate_ssm_context <- function(maf,
                                 motif = "WRCY",
                                 index = 3,
                                 genome_build,
                                 fastaPath,
                                 aid_alleles = FALSE,
                                 motif_strand = FALSE){

    if(missing(genome_build)) genome_build = NULL
    if(is.null(genome_build) && "maf_data" %in% class(maf)){
        genome_build = get_genome_build(maf)
    }
    index = rep_len(index, length(motif))
    for(i in seq_along(motif)){
        if(index[i] < 1 || index[i] > nchar(motif[i])){
            stop("index ", index[i], " is outside motif ", motif[i])
        }
    }
    needed = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
    if(!all(needed %in% colnames(maf))){
        stop("maf must contain the columns ", paste(needed, collapse = ", "))
    }
    fetch = .reference_fetcher(genome_build, if(missing(fastaPath)) NULL else fastaPath)

    ref = toupper(as.character(maf$Reference_Allele))
    alt = toupper(as.character(maf$Tumor_Seq_Allele2))
    is_snv = nchar(ref) == 1 & nchar(alt) == 1 & ref %in% c("A", "C", "G", "T") &
        alt %in% c("A", "C", "G", "T") & ref != alt

    # one reference lookup per distinct position, wide enough for every motif
    flank = max(nchar(motif)) - 1
    chrom = as.character(maf$Chromosome)
    pos = as.numeric(maf$Start_Position)
    key = paste(chrom, pos, sep = ":")
    ukey = unique(key[is_snv])
    win = rep(NA_character_, length(key))
    if(length(ukey)){
        uchrom = sub(":[^:]*$", "", ukey)
        upos = as.numeric(sub(".*:", "", ukey))
        useq = fetch(uchrom, upos - flank, upos + flank)
        win[is_snv] = useq[match(key[is_snv], ukey)]
    }
    centre = flank + 1

    for(i in seq_along(motif)){
        maf[[motif[i]]] = .motif_status(win, ref, alt, is_snv, motif[i], index[i], centre,
                                         aid_alleles, strand_col = FALSE)
        if(motif_strand){
            maf[[paste0(motif[i], "_strand")]] = .motif_status(win, ref, alt, is_snv, motif[i],
                                                                 index[i], centre, aid_alleles,
                                                                 strand_col = TRUE)
        }
    }
    l1 = nchar(motif[1])
    maf$seq = ifelse(is.na(win), NA_character_, substr(win, centre - l1 + 1, centre + l1 - 1))
    maf
}


#' @title Enrichment of mutations at a sequence motif
#'
#' @description Tests whether the single-nucleotide variants in each gene or
#' region fall at the target base of a motif (e.g. the C of the AID hotspot
#' WRCY) more often than expected from that region's own sequence.
#'
#' @details Only SNVs whose reference base could be a motif target are
#' counted (for WRCY: C on the forward strand or G on the reverse, i.e. any
#' C/G). The observed fraction is the share of those SNVs labelled SITE by
#' [annotate_ssm_context()]. The expected fraction is the share of the same
#' bases in the region's reference sequence that are motif targets, so the
#' comparison controls for base composition and for C/G mutations being more
#' common than A/T ones. A one-sided binomial test asks whether the observed
#' fraction exceeds the expected one.
#'
#' With `aid_alleles = TRUE`, SITE also requires an AID-type allele and the
#' expected fraction is multiplied by 2/3 (two of the three possible
#' substitutions), which assumes the three alternatives are equally likely
#' at random.
#'
#' The expected fraction comes from the whole region's sequence, so the SNVs
#' should come from data that cover it evenly: use genome data, or restrict
#' `regions` to the captured sequence for exome/capture data.
#'
#' Regions: supply `regions` (data frame with columns chrom, start, end, name;
#' mutations are assigned to every region they fall in). Otherwise each value
#' of `group_col` (default Hugo_Symbol) is one region spanning its SNVs, padded
#' by `pad` bp on each side; this assumes one chromosome per group.
#'
#' @param maf MAF data frame (columns as for [annotate_ssm_context()]).
#' @param motif Motif in IUPAC code (default "WRCY").
#' @param index Position of the target base in the motif (default 3).
#' @param genome_build "grch37" or "hg38" (taken from a maf_data object).
#' @param fastaPath Optional path to an indexed FASTA of the matching build.
#' @param regions Optional data frame of regions (chrom, start, end, name).
#' @param group_col Column defining groups when `regions` is not given.
#' @param pad Padding (bp) added to each group's span. Default 0.
#' @param aid_alleles Require AID-type alleles at SITE (default FALSE).
#' @param min_snv Groups with fewer eligible SNVs are reported without a test.
#'   Default 5.
#'
#' @return A data frame with one row per region: name, chrom, start, end,
#'   n_snv (eligible SNVs), n_site, observed_fraction, expected_fraction,
#'   enrichment (observed / expected) and p_value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' maf = GAMBLR.open::get_ssm_by_regions(regions_list = c("12:7055631-7070479", "18:60790579-60987361"))
#' calc_motif_enrichment(maf, genome_build = "grch37")
#' }
calc_motif_enrichment <- function(maf,
                                  motif = "WRCY",
                                  index = 3,
                                  genome_build,
                                  fastaPath,
                                  regions,
                                  group_col = "Hugo_Symbol",
                                  pad = 0,
                                  aid_alleles = FALSE,
                                  min_snv = 5){

    if(missing(genome_build)) genome_build = NULL
    if(is.null(genome_build) && "maf_data" %in% class(maf)){
        genome_build = get_genome_build(maf)
    }
    fp = if(missing(fastaPath)) NULL else fastaPath
    fetch = .reference_fetcher(genome_build, fp)
    ann = annotate_ssm_context(as.data.frame(maf), motif = motif, index = index,
                               genome_build = genome_build,
                               fastaPath = if(is.null(fp)) NULL else fp,
                               aid_alleles = aid_alleles)
    status = ann[[motif]]
    target = .motif_target_bases(motif, index)   # e.g. c("C", "G") for WRCY
    eligible = !is.na(status) & toupper(ann$Reference_Allele) %in% target
    chrom = sub("^chr", "", as.character(ann$Chromosome))
    pos = as.numeric(ann$Start_Position)

    if(missing(regions)){
        if(!group_col %in% colnames(ann)) stop("group_col ", group_col, " is not a column of maf")
        g = as.character(ann[[group_col]])
        keep = eligible & !is.na(g) & g != ""
        if(!any(keep)) stop("no eligible SNVs to test")
        regions = do.call(rbind, lapply(split(which(keep), g[keep]), function(ix){
            data.frame(chrom = names(sort(table(chrom[ix]), decreasing = TRUE))[1],
                       start = min(pos[ix]) - pad, end = max(pos[ix]) + pad,
                       name = g[ix[1]], stringsAsFactors = FALSE)
        }))
    }else{
        regions = as.data.frame(regions)
        if(!all(c("chrom", "start", "end", "name") %in% colnames(regions))){
            stop("regions needs columns chrom, start, end and name")
        }
        regions$chrom = sub("^chr", "", as.character(regions$chrom))
    }
    regions$start = pmax(1, regions$start)

    seqs = fetch(regions$chrom, regions$start, regions$end)
    out = lapply(seq_len(nrow(regions)), function(i){
        inr = eligible & chrom == regions$chrom[i] & pos >= regions$start[i] & pos <= regions$end[i]
        n_snv = sum(inr)
        n_site = sum(status[inr] == "SITE")
        exp_frac = .expected_site_fraction(seqs[i], motif, index)
        if(aid_alleles) exp_frac = exp_frac * 2 / 3
        p = if(n_snv >= min_snv && is.finite(exp_frac) && exp_frac > 0){
            stats::binom.test(n_site, n_snv, p = exp_frac, alternative = "greater")$p.value
        }else NA_real_
        data.frame(name = regions$name[i], chrom = regions$chrom[i], start = regions$start[i],
                   end = regions$end[i], n_snv = n_snv, n_site = n_site,
                   observed_fraction = if(n_snv) n_site / n_snv else NA_real_,
                   expected_fraction = exp_frac,
                   enrichment = if(n_snv && exp_frac > 0) (n_site / n_snv) / exp_frac else NA_real_,
                   p_value = p, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}


# ---- internal helpers -------------------------------------------------------

.iupac_regex <- c(A = "A", C = "C", G = "G", T = "T", R = "[AG]", Y = "[CT]", S = "[GC]",
                  W = "[AT]", K = "[GT]", M = "[AC]", B = "[CGT]", D = "[AGT]", H = "[ACT]",
                  V = "[ACG]", N = "[ACGT]")
.iupac_comp <- c(A = "T", T = "A", C = "G", G = "C", R = "Y", Y = "R", S = "S", W = "W",
                 K = "M", M = "K", B = "V", D = "H", H = "D", V = "B", N = "N")

.motif_regex <- function(m){
    b = strsplit(toupper(m), "")[[1]]
    if(!all(b %in% names(.iupac_regex))) stop("motif ", m, " contains non-IUPAC letters")
    paste0("^", paste(.iupac_regex[b], collapse = ""), "$")
}
.revcomp_motif <- function(m) paste(rev(.iupac_comp[strsplit(toupper(m), "")[[1]]]), collapse = "")

# reference bases that can be the target: the forward target letter and its complement
.motif_target_bases <- function(motif, index){
    b = substr(toupper(motif), index, index)
    fwd = strsplit(gsub("[][]", "", .iupac_regex[b]), "")[[1]]
    rev = strsplit(gsub("[][]", "", .iupac_regex[.iupac_comp[b]]), "")[[1]]
    unique(c(fwd, rev))
}

# SITE/MOTIF/FALSE (or the strand of SITE) for each row
.motif_status <- function(win, ref, alt, is_snv, motif, index, centre, aid_alleles, strand_col){
    L = nchar(motif)
    rc = .revcomp_motif(motif)
    rindex = L - index + 1
    fwd_re = .motif_regex(motif)
    rev_re = .motif_regex(rc)
    ok = is_snv & !is.na(win)
    fwd_sub = substr(win, centre - index + 1, centre - index + L)
    rev_sub = substr(win, centre - rindex + 1, centre - rindex + L)
    fwd_site = ok & grepl(fwd_re, fwd_sub)
    rev_site = ok & grepl(rev_re, rev_sub)
    if(aid_alleles){
        comp = c(A = "T", C = "G", G = "C", T = "A")
        allowed = function(r, a) (r == "C" & a %in% c("T", "G")) | (r == "A" & a %in% c("G", "T"))
        fwd_site = fwd_site & allowed(ref, alt)
        rev_site = rev_site & allowed(unname(comp[ref]), unname(comp[alt]))
    }
    if(strand_col){
        return(ifelse(fwd_site & rev_site, "both", ifelse(fwd_site, "+", ifelse(rev_site, "-", NA_character_))))
    }
    in_motif = rep(FALSE, length(win))
    for(s in seq_len(L)){   # every L-mer that contains the centre base
        sub = substr(win, centre - L + s, centre - 1 + s)
        in_motif = in_motif | (ok & (grepl(fwd_re, sub) | grepl(rev_re, sub)))
    }
    out = ifelse(fwd_site | rev_site, "SITE", ifelse(in_motif, "MOTIF", "FALSE"))
    out[!ok] = NA_character_
    out
}

# fraction of target-eligible bases in seq that are motif targets (either strand)
.expected_site_fraction <- function(seq, motif, index){
    if(is.na(seq) || !nchar(seq)) return(NA_real_)
    s = Biostrings::DNAString(toupper(seq))
    L = nchar(motif)
    f = Biostrings::matchPattern(Biostrings::DNAString(motif), s, fixed = FALSE)
    r = Biostrings::matchPattern(Biostrings::DNAString(.revcomp_motif(motif)), s, fixed = FALSE)
    targets = unique(c(IRanges::start(f) + index - 1, IRanges::start(r) + (L - index)))
    bases = strsplit(as.character(s), "")[[1]]
    n_eligible = sum(bases %in% .motif_target_bases(motif, index))
    if(n_eligible == 0) return(NA_real_)
    length(targets) / n_eligible
}

# returns function(chrom, start, end) -> character sequences, matching the
# reference's chromosome naming whatever the input uses
.reference_fetcher <- function(genome_build, fastaPath = NULL){
    if(!is.null(fastaPath)){
        if(!file.exists(fastaPath)) stop("fastaPath ", fastaPath, " does not exist")
        fa = Rsamtools::FaFile(fastaPath)
        seqnames = GenomeInfoDb::seqlevels(Rsamtools::scanFaIndex(fa))
        getter = function(gr) as.character(Rsamtools::getSeq(fa, gr))
    }else{
        if(is.null(genome_build)){
            stop("Provide genome_build (or a maf_data object that carries it) or fastaPath")
        }
        pkg = switch(genome_build, grch37 = "BSgenome.Hsapiens.UCSC.hg19",
                     hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                     stop("unsupported genome_build: ", genome_build))
        if(!pkg %in% BSgenome::installed.genomes()){
            stop("No fastaPath given and ", pkg, " is not installed. Install it or supply fastaPath.")
        }
        genome = BSgenome::getBSgenome(pkg)
        seqnames = GenomeInfoDb::seqnames(genome)
        getter = function(gr) as.character(Biostrings::getSeq(genome, gr))
    }
    use_chr = any(grepl("^chr", seqnames))
    function(chrom, start, end){
        chrom = sub("^chr", "", as.character(chrom))
        if(use_chr) chrom = paste0("chr", chrom)
        out = rep(NA_character_, length(chrom))
        known = chrom %in% seqnames
        if(any(known)){
            gr = GenomicRanges::GRanges(chrom[known], IRanges::IRanges(pmax(1, start[known]), end[known]))
            out[known] = getter(gr)
        }
        out
    }
}
