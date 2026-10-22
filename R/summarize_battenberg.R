#' @title Summarize allele-specific copy number from Battenberg
#'
#' @description Summarize Battenberg allele-specific copy number calls for
#' one or many samples: purity, ploidy and fit quality; whole-genome doubling;
#' loss of heterozygosity (LOH) split by mechanism (copy-neutral, with copy
#' loss, with copy gain); allelic imbalance without LOH; genomic instability;
#' and LOH at recurrent loci. Returns a per-sample summary and a table of the
#' LOH segments.
#'
#' @details
#' \strong{Inputs.} For each sample, `<prefix>_subclones.txt` is read along
#' with, when present in the same directory, `<prefix>_cellularity_ploidy.txt`,
#' `<prefix>_rho_and_psi.txt` and `<prefix>_copynumber_solution_status.txt`.
#' The sample name is the file prefix. Only autosomes are used (sex
#' chromosomes depend on sex). Chromosome names are compared without a "chr"
#' prefix.
#'
#' \strong{Copy number states.} Each segment's state is the clonal state
#' (`nMaj1_A`, `nMin1_A`) or, for subclonal segments, whichever of the two
#' states covers more cells. Segments with `frac1_A < 1` count as subclonal.
#' The sample's baseline is its modal total copy number (by length).
#'
#' \strong{LOH classes} (minor allele = 0, total > 0), relative to the
#' baseline: `loh_neutral` (total = baseline; e.g. 2+0 in a diploid genome,
#' 4+0 in a tetraploid one), `loh_loss` (total below baseline; e.g. 1+0 in a
#' diploid genome, 2+0 in a tetraploid one) and `loh_gain` (total above
#' baseline). Total copy number 0 is a homozygous deletion. Segments of the
#' same class are merged into events, including across interruptions of up to
#' `merge_max_gap` (e.g. a short segment of another state inside a large LOH
#' region), before counting events and measuring their sizes; genome fractions
#' use all segments.
#'
#' \strong{Terminal events} reach the first or last called position of their
#' chromosome; the others are interstitial (counted in
#' `loh_neutral_interstitial_events` and `loh_loss_interstitial_events`, for
#' events of at least `min_event_size`). Copy-neutral LOH from mitotic recombination is typically
#' terminal, whereas LOH from deletion is often interstitial.
#'
#' \strong{Whole-genome doubling} (`doubled`) is a heuristic, not an
#' evolutionary inference: TRUE when at least half of the autosomal genome has
#' a major allele copy number of 2 or more. It is distinct from the average
#' ploidy (`ploidy`, from Battenberg; `mean_total_cn`, from the calls) and from
#' the modal baseline (`modal_cn`, `baseline_cn`) used to call relative gains
#' and losses. In doubled
#' genomes, LOH with an even major allele copy number (2+0, 4+0) is
#' consistent with LOH before doubling (a single-copy or copy-neutral LOH that
#' was then doubled), and an odd one (1+0, 3+0) with loss after doubling.
#'
#' \strong{Shares.} `loh_neutral_share` is the copy-neutral fraction of LOH,
#' `loh_neutral_fraction / (loh_neutral_fraction + loh_loss_fraction)`.
#' `loss_fraction` and `gain_fraction` are the genome below and above the
#' baseline (including homozygous deletions in the former), and
#' `loss_with_loh_share` is the share of that lost genome with LOH,
#' `loh_loss_fraction / loss_fraction`. All use the same calls and baseline;
#' they are NA when the denominator is 0.
#'
#' \strong{Event counts and small events.} Event counts and sizes include
#' only merged events of at least `min_event_size`; genome fractions include
#' every segment. A sample can therefore have a small nonzero LOH fraction
#' with zero counted events (`loh_neutral_small_events` and
#' `loh_loss_small_events` count the events below the threshold).
#'
#' \strong{Allele-specific properties.} Three properties are recorded for
#' each segment independently, in `segments`: `cn_change` (total copy number
#' below, at or above the baseline: loss, neutral, gain), `homologues` (both
#' parental homologues present, i.e. minor allele > 0; one, i.e. complete
#' LOH; none, i.e. total 0) and `balance` (for segments with both homologues:
#' balanced, major = minor, or imbalanced). For example, at a baseline of 4,
#' 1+1 is a loss retaining both homologues (balanced), 2+0 a loss with
#' complete LOH, 2+1 a loss retaining both (imbalanced), 3+1 copy-neutral
#' imbalance and 4+0 copy-neutral complete LOH. The allele copy numbers and
#' the `allelic_state` label (e.g. "2+1") are kept. Allele copy numbers must
#' be integers (as Battenberg reports them); other values are unknown, not
#' rounded. For a subclonal segment, each property comes from the dominant
#' state when it covers at least `min_clonal_fraction` of tumour cells,
#' otherwise only if both states agree; when they disagree the property is
#' NA. The baseline for these properties, `baseline_cn`, is the modal total
#' copy number over segments whose total is known; it normally equals
#' `modal_cn` (which uses each segment's majority state).
#'
#' \strong{Allele-specific fractions.} Their denominator is the assessable
#' genome: autosomal segments whose copy number change, homologues and (for
#' segments with both homologues) balance are all known
#' (`assessable_genome_mb`; `assessable_fraction` is its share of the called
#' genome, `called_genome_mb`, and `coverage_fraction` the called share of
#' the reference autosomes, `reference_autosome_mb`, excluding
#' centromeres). Unassessable or uncalled regions are not counted as
#' heterozygous or copy-neutral, and these fractions are NA when nothing is
#' assessable. Fractions of the assessable genome: `complete_loh_fraction`
#' (one homologue) and its split by copy number change,
#' `complete_loh_loss/neutral/gain_fraction`; `relative_loss_fraction`;
#' `neutral_fraction`; `no_homologue_fraction` (total 0).
#' `loh_copy_neutral_share` is the copy-neutral share of complete LOH
#' (including LOH with gain). `allelic_states` gives the genome in each
#' allele-specific state per sample (e.g. 2+2 vs 4+0), from which other
#' combinations can be derived; it covers segments whose exact state is known,
#' so its fractions (of the assessable genome) can sum to slightly less than 1.
#' `mean_total_cn` is the length-weighted mean total copy number over segments
#' whose exact total is known.
#'
#' These measure detectable changes in parental allele composition only. A
#' repair or recombination event that copies one homologue from another copy
#' of the same parental homologue (e.g. a sister chromatid, or after
#' doubling, the other copy of that parent: AA'BB' to AA'B'B') leaves the
#' allele state unchanged (2+2) and is invisible here. They are endpoint
#' states, not rates of loss or recombination.
#'
#' \strong{Fitting stages and QC.} From `<prefix>_rho_and_psi.txt`, by row:
#' `ascat_purity`/`ascat_ploidy` (initial ASCAT fit), `clonal_fit_purity`,
#' `clonal_fit_ploidy` and `clonal_fit_score` (the FRAC_GENOME clonal refit,
#' which the final profile is fitted from, and its distance) and
#' `ref_seg_ploidy`; `final_profile_ploidy` is the ploidy in
#' `<prefix>_cellularity_ploidy.txt` (the same as `ploidy`). The
#' `ploidy_delta_*` columns are signed differences (later stage minus earlier).
#' `fit_distance` is unchanged: the distance of the row marked `is.best`.
#' `refit_suggestions` counts the rows of `<prefix>_refit_suggestion.txt` (NA
#' if absent). Missing or malformed files leave these NA. QC annotation
#' (`allelic_qc_status`, `allelic_qc_reasons` and `qc_*` flags) is added by
#' [battenberg_qc].
#'
#' \strong{Homozygous deletion size.} `homdel_mb` is the total size of
#' homozygous deletion (total copy number 0) and `largest_homdel_mb` the
#' largest run of consecutive homozygously deleted segments (deleted sequence
#' only; unlike `loh_segments` events, not merged across gaps).
#' Homozygous deletions of tens of Mb are not compatible with a viable cell
#' and usually point to a misfit solution; [battenberg_qc] flags them.
#'
#' \strong{Excluded LOH.} `excluded_loh_fraction` is the share of the called
#' genome that is LOH by majority state but not assessable (e.g. subclonal
#' segments whose states disagree). A low `complete_loh_fraction` with a high
#' `excluded_loh_fraction` means LOH was excluded, not absent; the original
#' `loh_*` columns report LOH by majority state over the called genome.
#'
#' \strong{Complete LOH events} (`complete_loh_events`, and the `cnloh_*`
#' summary columns). Built from the allele-specific properties (assessable
#' segments with one homologue and an exact total copy number), so
#' unresolved subclonal segments are left out; consecutive segments with the
#' same total copy number form an event. Unlike `loh_segments` (merged by
#' majority state over the called genome), these are allele-specific.
#' Each event's `loh_type` (neutral, loss or gain) is relative to
#' `cnloh_reference`. What counts as copy-neutral is
#' set by `cnloh_reference`:
#' \itemize{
#'   \item "sample" (default): total copy number equals the sample's baseline
#'     (`baseline_cn`), e.g. 2+0 in a diploid genome and 4+0 in a doubled one.
#'   \item "diploid": total copy number 2 whatever the baseline, e.g. 2+0 in a
#'     doubled genome too (where 4+0 is then LOH with gain).
#' }
#' In a doubled genome, 4+0 is most simply copy-neutral LOH before doubling
#' (2+0, then doubled), while 2+0 is most simply a single-copy loss before
#' doubling (1+0, then doubled); a single copy-neutral event after doubling
#' gives 3+1, not LOH. The event table holds every complete LOH event, with
#' `neutral_vs_sample` and `neutral_vs_diploid` so either definition can be
#' applied without re-running; `large` marks events of at least
#' `min_event_size` and `counted` the large copy-neutral ones under the
#' chosen reference (those in the `cnloh_*` counts).
#' Events are merged across uncalled sequence (e.g. centromeres) and across
#' up to `merge_max_gap` of segments in other states.
#' Each event's `extent`, measured against the called part of its chromosome
#' (uncalled regions, e.g. near centromeres, don't count against it), is
#' `whole_chromosome` (at least `arm_fraction` of the chromosome),
#' `arm` (at least `arm_fraction` of an arm), `terminal` (reaching the first
#' or last called position) or `interstitial`. These extents suggest
#' different mechanisms (e.g. missegregation with reduplication for whole
#' chromosomes, mitotic recombination for terminal events), so whole-
#' chromosome events are reported separately: `cnloh_fraction` includes them
#' and `cnloh_segmental_fraction` excludes them (both of the assessable
#' genome, from all segments); counts and Mb per extent use the counted
#' events. Two copy-neutral shares of LOH leave out whole-chromosome LOH of
#' any copy number (e.g. monosomy as well as whole-chromosome copy-neutral
#' LOH), from both numerator and denominator:
#' `loh_copy_neutral_share_segmental` (by length: copy-neutral share of
#' complete LOH with an exact total copy number) and
#' `loh_copy_neutral_event_share` (by events of at least `min_event_size`,
#' each event counting once, so not dominated by a few large events); NA when
#' there is no such LOH. Both follow `cnloh_reference`, unlike
#' `loh_copy_neutral_share`, which is always relative to the sample's
#' baseline and includes whole-chromosome events. Also per event: `subclonal_fraction` (share in subclonal
#' segments), `baf_deviation` (mean |observed - expected| major-allele BAF at
#' the fitted purity; large values suggest a poor fit) and the `loh_genes` it
#' covers. The older `loh_neutral_*` columns use majority states over the
#' called genome and the sample baseline.
#'
#' \strong{Allelic imbalance without LOH} (minor allele > 0, major > minor) is
#' split the same way: `imbalance_gain` (e.g. 2+1 in a diploid genome),
#' `imbalance_neutral` (e.g. 3+1 in a tetraploid one) and `imbalance_loss`.
#'
#' @param battenberg_files Optional table of samples and file paths, e.g. from
#' `GAMBLR.results::get_battenberg()`: columns sample_id, subclones_file and
#' projection (the genome build of each sample), and optionally
#' cellularity_ploidy_file, rho_and_psi_file and solution_status_file.
#' Samples whose subclones file is missing are skipped with a message. When
#' given, `battenberg_dir`, `subclones_files` and `genome_build` are not used.
#' @param battenberg_dir Directory to search (recursively) for
#' `*_subclones.txt` files; the sample name is the file prefix and the other
#' files are looked up next to it. Ignored when `subclones_files` is given.
#' @param subclones_files Optional vector of `*_subclones.txt` paths.
#' @param genome_build "grch37" (default) or "hg38"; used for chromosome arm
#' coordinates and the LOH loci (with `battenberg_dir` or `subclones_files`).
#' @param min_event_size Minimum size (bp) of a merged event for it to be
#' counted in the event counts and size statistics. Default 1e6. Genome
#' fractions use all segments.
#' @param merge_max_gap Largest gap (bp) between two events of the same class
#' on a chromosome that is merged over. Default 2e6.
#' @param arm_fraction Fraction of a chromosome arm an LOH event must cover to
#' count as arm-level LOH. Default 0.9.
#' @param min_clonal_fraction For subclonal segments, the fraction of tumour
#' cells the dominant state must cover for its allele-specific properties to
#' be used when the two states disagree. Default 0.8.
#' @param cnloh_reference What copy-neutral LOH means for the `cnloh_*`
#' summary columns and the `counted` events: "sample" (default; total equal
#' to the sample's baseline copy number) or "diploid" (total copy number 2).
#' See Details.
#' @param qc_thresholds,clonality_dist_metric Passed to [battenberg_qc] to
#' annotate the samples with allele-specific QC (see there). Defaults: the
#' provisional thresholds and Battenberg's default distance metric (0).
#' @param loh_genes Genes whose LOH status is reported (one column each).
#' Default HLA-A, HLA-B, B2M, TP53, TNFRSF14 and CDKN2A. Set to NULL to skip.
#' @param verbose Print progress. Default FALSE.
#'
#' @return A list with \code{summary} (one row per sample; see Details for
#' the metric definitions), \code{loh_segments} (one row per merged LOH or
#' homozygous deletion event, with its class, size, arm and whether it is
#' terminal) and \code{segments} (every autosomal Battenberg segment with its
#' observed BAF, LogR and LogR-based total copy number `ntot`, its clonal
#' `major`/`minor`/`total` copy number, whether it is subclonal, and its
#' `state`: loh_neutral, loh_loss, loh_gain, homdel, imbalance_* or balanced).
#' The `segments` table is what [GAMBLR.viz::plot_battenberg] plots; it also
#' has the allele-specific properties (`major_cn`, `minor_cn`, `total_cn`,
#' `allelic_state`, `cn_change`, `homologues`, `balance`, `assessable`).
#' Also \code{allelic_states} (genome per allele-specific state per sample)
#' and \code{complete_loh_events} (allele-specific complete LOH events with
#' their extent, allele state and type under `cnloh_reference`), e.g. for
#' [battenberg_to_loh_status]. See Details.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # all GAMBL genomes, via the GAMBLR config
#' bb_files = GAMBLR.results::get_battenberg(GAMBLR.results::get_gambl_metadata())
#' bb = summarize_battenberg(bb_files)
#'
#' # or every result under a directory
#' bb = summarize_battenberg(battenberg_dir = "/path/to/battenberg/results")
#' bb$summary
#' write.table(bb$summary, "battenberg_summary.tsv", sep = "\t",
#'             quote = FALSE, row.names = FALSE)
#' }
summarize_battenberg = function(battenberg_files = NULL,
                                battenberg_dir = NULL,
                                subclones_files = NULL,
                                genome_build = "grch37",
                                min_event_size = 1e6,
                                merge_max_gap = 2e6,
                                arm_fraction = 0.9,
                                min_clonal_fraction = 0.8,
                                cnloh_reference = "sample",
                                qc_thresholds = NULL,
                                clonality_dist_metric = 0,
                                loh_genes = c("HLA-A", "HLA-B", "B2M", "TP53",
                                              "TNFRSF14", "CDKN2A"),
                                verbose = FALSE){
  cnloh_reference = match.arg(cnloh_reference, c("sample", "diploid"))
  # one row per sample: its name, build and files
  if(!is.null(battenberg_files)){
    if(is.character(battenberg_files) && is.null(dim(battenberg_files))){
      # a directory given as the first argument
      battenberg_dir = battenberg_files
      battenberg_files = NULL
    }
  }
  if(!is.null(battenberg_files)){
    jobs = as.data.frame(battenberg_files)
    needed = c("sample_id", "subclones_file", "projection")
    missing_cols = setdiff(needed, colnames(jobs))
    if(length(missing_cols) > 0){
      stop(paste("battenberg_files is missing columns:",
                 paste(missing_cols, collapse = ", "),
                 "(use GAMBLR.results::get_battenberg())"))
    }
    jobs$sample = jobs$sample_id
    jobs$build = jobs$projection
    absent = !file.exists(jobs$subclones_file)
    if(any(absent)){
      message(sum(absent), " sample(s) have no subclones file and are skipped")
      jobs = jobs[!absent, , drop = FALSE]
    }
  }else{
    if(is.null(subclones_files)){
      if(is.null(battenberg_dir)){
        stop("give battenberg_files, battenberg_dir or subclones_files")
      }
      subclones_files = list.files(battenberg_dir, pattern = "_subclones\\.txt$",
                                   recursive = TRUE, full.names = TRUE)
    }
    sample = sub("_subclones\\.txt$", "", basename(subclones_files))
    prefix = file.path(dirname(subclones_files), sample)
    jobs = data.frame(sample = sample, subclones_file = subclones_files,
                      cellularity_ploidy_file = paste0(prefix, "_cellularity_ploidy.txt"),
                      rho_and_psi_file = paste0(prefix, "_rho_and_psi.txt"),
                      solution_status_file = paste0(prefix, "_copynumber_solution_status.txt"),
                      build = genome_build, stringsAsFactors = FALSE)
  }
  if(nrow(jobs) == 0){
    stop("no *_subclones.txt files found")
  }
  jobs$build = ifelse(tolower(jobs$build) %in% c("grch37", "hg19", "hs37d5"), "grch37",
                      ifelse(tolower(jobs$build) %in% c("hg38", "grch38"), "hg38", NA))
  if(anyNA(jobs$build)){
    stop("genome builds must be grch37 or hg38")
  }

  # chromosome arms and LOH loci for each build in use
  reference = lapply(stats::setNames(unique(jobs$build), unique(jobs$build)), function(build){
    arms = as.data.frame(if(build == "grch37") GAMBLR.data::chromosome_arms_grch37
                         else GAMBLR.data::chromosome_arms_hg38)
    arms$chromosome = gsub("^chr", "", as.character(arms$chromosome))
    loci = NULL
    if(length(loh_genes) > 0){
      loci = as.data.frame(suppressMessages(
        gene_to_region(gene_symbol = loh_genes, projection = build, return_as = "df")))
      gene_col = grep("hugo|gene", colnames(loci), ignore.case = TRUE)[1]
      colnames(loci)[seq_len(3)] = c("chrom", "start", "end")
      loci$chrom = gsub("^chr", "", as.character(loci$chrom))
      loci$gene = loh_genes[match(toupper(loci[[gene_col]]), toupper(loh_genes))]
    }
    list(arms = arms, loci = loci)
  })

  results = lapply(seq_len(nrow(jobs)), function(i){
    job = jobs[i, ]
    if(verbose){
      message("summarizing ", job$sample)
    }
    companion = function(column){
      if(column %in% colnames(job)) job[[column]] else NA_character_
    }
    tryCatch(.summarize_one_battenberg(
      job$subclones_file, sample = job$sample,
      cellularity_file = companion("cellularity_ploidy_file"),
      rho_file = companion("rho_and_psi_file"),
      status_file = companion("solution_status_file"),
      arms = reference[[job$build]]$arms, loci = reference[[job$build]]$loci,
      min_event_size = min_event_size, merge_max_gap = merge_max_gap,
      arm_fraction = arm_fraction, min_clonal_fraction = min_clonal_fraction,
      cnloh_reference = cnloh_reference),
      error = function(e){
        warning("skipping ", job$sample, ": ", conditionMessage(e))
        NULL
      })
  })
  results = results[!vapply(results, is.null, logical(1))]
  out = list(summary = do.call(rbind, lapply(results, `[[`, "summary")),
       loh_segments = do.call(rbind, lapply(results, `[[`, "loh_segments")),
       segments = do.call(rbind, lapply(results, `[[`, "segments")),
       allelic_states = do.call(rbind, lapply(results, `[[`, "allelic_states")),
       complete_loh_events = do.call(rbind, lapply(results, `[[`, "complete_loh_events")))
  battenberg_qc(out, qc_thresholds = qc_thresholds,
                clonality_dist_metric = clonality_dist_metric)
}

# One sample: returns list(summary = one-row data frame, loh_segments)
.summarize_one_battenberg = function(subclones_file, sample, cellularity_file,
                                     rho_file, status_file, arms, loci,
                                     min_event_size, merge_max_gap,
                                     arm_fraction, min_clonal_fraction,
                                     cnloh_reference = "sample"){
  x = utils::read.delim(subclones_file, stringsAsFactors = FALSE)
  x$chr = gsub("^chr", "", as.character(x$chr))
  x = x[x$chr %in% as.character(1:22), , drop = FALSE]
  x = x[order(as.integer(x$chr), x$startpos), , drop = FALSE]
  if(nrow(x) == 0){
    stop("no autosomal segments")
  }
  same_chrom = x$chr[-1] == x$chr[-nrow(x)]
  if(any(same_chrom & x$startpos[-1] < x$endpos[-nrow(x)])){
    stop("overlapping segments; lengths would be double counted")
  }
  len = as.numeric(x$endpos - x$startpos)

  # clonal (majority) state of each segment
  use2 = !is.na(x$frac2_A) & !is.na(x$frac1_A) & x$frac2_A > x$frac1_A
  maj = ifelse(use2, x$nMaj2_A, x$nMaj1_A)
  min = ifelse(use2, x$nMin2_A, x$nMin1_A)
  tot = maj + min
  called = !is.na(tot)
  genome = sum(len[called])
  frac = function(cond) sum(len[cond & called], na.rm = TRUE) / genome
  modal = as.numeric(names(which.max(tapply(len[called], tot[called], sum))))

  state = ifelse(!called, NA_character_,
           ifelse(tot == 0, "homdel",
           ifelse(min == 0, ifelse(tot == modal, "loh_neutral",
                              ifelse(tot < modal, "loh_loss", "loh_gain")),
           ifelse(maj != min, ifelse(tot == modal, "imbalance_neutral",
                                ifelse(tot < modal, "imbalance_loss",
                                       "imbalance_gain")),
                  "balanced"))))
  doubled = frac(maj >= 2) >= 0.5

  # allele-specific properties (Phase 1): copy number change relative to the
  # baseline, which parental homologues are present, and whether the retained
  # homologues are balanced; unknown (NA) when the subclonal states disagree
  alleles = .battenberg_allele_properties(x, min_clonal_fraction)

  # merged events: consecutive segments of the same LOH/homdel class
  first_pos = tapply(x$startpos, x$chr, min)
  last_pos = tapply(x$endpos, x$chr, max)
  event_class = ifelse(state %in% c("loh_neutral", "loh_loss", "loh_gain",
                                    "homdel"), state, NA_character_)
  events = .battenberg_events(x, event_class, sample, major = maj, total = tot,
                              merge_max_gap = merge_max_gap, arms = arms,
                              arm_fraction = arm_fraction,
                              first_pos = first_pos, last_pos = last_pos)
  sized = if(is.null(events)) events else events[events$size >= min_event_size, ]
  event_stats = function(cls, ev = events){
    sized = if(is.null(ev)) ev else ev[ev$size >= min_event_size, ]
    e = if(is.null(sized)) NULL else sized[sized$class == cls, ]
    n = if(is.null(e)) 0L else nrow(e)
    c(n = n,
      small = if(is.null(ev)) 0L else sum(ev$class == cls & ev$size < min_event_size),
      median_size_mb = if(n > 0) stats::median(e$size) / 1e6 else NA_real_,
      terminal_fraction = if(n > 0) mean(e$terminal) else NA_real_,
      interstitial = if(n > 0) sum(!e$terminal) else 0,
      arm_level = if(n > 0) sum(e$arm_level) else 0)
  }
  neutral = event_stats("loh_neutral")
  loss = event_stats("loh_loss")

  # chromosomes entirely LOH (any class), e.g. whole-chromosome UPD
  loh_any = state %in% c("loh_neutral", "loh_loss", "loh_gain")
  chrom_loh = tapply(len * loh_any, x$chr, sum) / tapply(len, x$chr, sum)

  # fit and solution quality, when the files are present
  # NULL when absent or unreadable (malformed), so fitting fields become NA
  read_if = function(f){
    if(is.na(f) || !file.exists(f)) return(NULL)
    tryCatch(utils::read.delim(f, stringsAsFactors = FALSE), error = function(e) NULL)
  }
  cp = read_if(cellularity_file)
  rp = read_if(rho_file)
  if(!is.null(rp) && !all(c("ASCAT", "FRAC_GENOME") %in% rownames(rp))){
    # unexpected layout: keep the fitting-stage fields NA rather than guess
    rp = NULL
  }
  # fit_distance: the distance of the solution marked is.best (unchanged)
  fit_distance = NA_real_
  if(!is.null(rp) && "is.best" %in% colnames(rp)){
    best = rp[as.character(rp$is.best) %in% c("TRUE", "True"), , drop = FALSE]
    if(nrow(best) > 0) fit_distance = as.numeric(best$distance[1])
  }
  # each fitting stage by name (rows ASCAT, FRAC_GENOME, REF_SEG); the final
  # profile is fitted from the FRAC_GENOME solution
  stage = function(row, column){
    if(is.null(rp) || !row %in% rownames(rp) || !column %in% colnames(rp)) return(NA_real_)
    suppressWarnings(as.numeric(rp[row, column]))
  }
  refit_file = sub("_subclones\\.txt$", "_refit_suggestion.txt", subclones_file)
  refit_suggestions = if(file.exists(refit_file)){
    r = tryCatch(utils::read.delim(refit_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if(is.null(r)) NA_integer_ else nrow(r)
  }else NA_integer_

  summary = data.frame(
    sample = sample,
    purity = if(!is.null(cp)) cp$cellularity[1] else NA_real_,
    ploidy = if(!is.null(cp)) cp$ploidy[1] else NA_real_,
    fit_distance = fit_distance,
    solution_status = if(!is.na(status_file) && file.exists(status_file))
      paste(readLines(status_file, warn = FALSE), collapse = " ") else NA_character_,
    modal_cn = modal,
    major_ge2_fraction = frac(maj >= 2),
    doubled = doubled,
    subclonal_fraction = frac(!is.na(x$frac1_A) & x$frac1_A < 1),
    loh_fraction = frac(loh_any),
    loh_neutral_fraction = frac(state == "loh_neutral"),
    loh_loss_fraction = frac(state == "loh_loss"),
    loh_gain_fraction = frac(state == "loh_gain"),
    # share of LOH that is copy-neutral rather than through copy loss
    loh_neutral_share = {
      n = frac(state == "loh_neutral"); l = frac(state == "loh_loss")
      if(n + l > 0) n / (n + l) else NA_real_
    },
    # genome below / above the baseline, and the share of lost genome with LOH
    loss_fraction = frac(tot < modal),
    gain_fraction = frac(tot > modal),
    loss_with_loh_share = {
      l = frac(tot < modal)
      if(l > 0) frac(state == "loh_loss") / l else NA_real_
    },
    loh_neutral_events = neutral[["n"]],
    loh_neutral_median_size_mb = neutral[["median_size_mb"]],
    loh_neutral_terminal_fraction = neutral[["terminal_fraction"]],
    loh_neutral_interstitial_events = neutral[["interstitial"]],
    loh_neutral_arm_level = neutral[["arm_level"]],
    loh_loss_events = loss[["n"]],
    loh_loss_median_size_mb = loss[["median_size_mb"]],
    loh_loss_terminal_fraction = loss[["terminal_fraction"]],
    loh_loss_interstitial_events = loss[["interstitial"]],
    loh_loss_arm_level = loss[["arm_level"]],
    whole_chromosome_loh = sum(chrom_loh >= 0.95, na.rm = TRUE),
    # LOH timing in doubled genomes (NA otherwise)
    loh_even_major_fraction = if(doubled) frac(min == 0 & tot > 0 & maj %% 2 == 0) else NA_real_,
    loh_odd_major_fraction = if(doubled) frac(min == 0 & tot > 0 & maj %% 2 == 1) else NA_real_,
    imbalance_gain_fraction = frac(state == "imbalance_gain"),
    imbalance_neutral_fraction = frac(state == "imbalance_neutral"),
    imbalance_loss_fraction = frac(state == "imbalance_loss"),
    homdel_events = if(is.null(events)) 0L else sum(events$class == "homdel"),
    homdel_fraction = frac(state == "homdel"),
    segments_per_gb = sum(called) / (genome / 1e9),
    stringsAsFactors = FALSE
  )

  # LOH class at each locus (NA when the locus isn't covered)
  if(!is.null(loci)){
    for(i in seq_len(nrow(loci))){
      hit = which(x$chr == loci$chrom[i] & x$startpos <= loci$end[i] &
                    x$endpos >= loci$start[i])
      value = if(length(hit) == 0) NA_character_ else {
        s = state[hit]
        loh = s[s %in% c("loh_neutral", "loh_loss", "loh_gain", "homdel")]
        if(length(loh) > 0) loh[1] else "none"
      }
      summary[[paste0("loh_", loci$gene[i])]] = value
    }
  }
  # allele-specific summary on the assessable genome
  a = alleles
  assessable = !is.na(a$cn_change) & !is.na(a$homologues) &
    (a$homologues != "both" | !is.na(a$balance))
  A = sum(len[assessable])
  afrac = function(cond) if(A > 0) sum(len[assessable & cond], na.rm = TRUE) / A else NA_real_
  share = function(num, den) if(!is.na(den) && den > 0) num / den else NA_real_
  is_ = function(v, value) !is.na(v) & v == value
  is_loss = is_(a$cn_change, "loss"); is_neutral = is_(a$cn_change, "neutral")
  is_gain = is_(a$cn_change, "gain")
  one = is_(a$homologues, "one"); none = is_(a$homologues, "none")
  autosome_mb = sum(arms$end[arms$chromosome %in% as.character(1:22)] -
                      arms$start[arms$chromosome %in% as.character(1:22)]) / 1e6
  totals = a$total_cn[assessable]
  final_ploidy = if(!is.null(cp)) cp$ploidy[1] else NA_real_
  ascat_ploidy = stage("ASCAT", "ploidy")
  clonal_ploidy = stage("FRAC_GENOME", "ploidy")
  extra = data.frame(
    ascat_purity = stage("ASCAT", "rho"),
    ascat_ploidy = ascat_ploidy,
    clonal_fit_purity = stage("FRAC_GENOME", "rho"),
    clonal_fit_ploidy = clonal_ploidy,
    clonal_fit_score = stage("FRAC_GENOME", "distance"),
    ref_seg_ploidy = stage("REF_SEG", "ploidy"),
    final_profile_ploidy = final_ploidy,
    ploidy_delta_ascat_clonal = clonal_ploidy - ascat_ploidy,
    ploidy_delta_clonal_final = final_ploidy - clonal_ploidy,
    ploidy_delta_ascat_final = final_ploidy - ascat_ploidy,
    refit_suggestions = refit_suggestions,
    # LOH (by majority state) in segments excluded as unassessable, as a
    # share of the called genome
    excluded_loh_fraction = sum(len[state %in% c("loh_neutral", "loh_loss", "loh_gain") &
                                      !assessable], na.rm = TRUE) / genome,
    baseline_cn = a$baseline,
    # length-weighted mean over segments whose exact total is known (a
    # subclonal segment can be assessable without an exact total)
    mean_total_cn = {
      known_total = !is.na(a$total_cn)
      if(any(known_total)) sum(a$total_cn[known_total] * len[known_total]) / sum(len[known_total]) else NA_real_
    },
    reference_autosome_mb = autosome_mb,
    called_genome_mb = genome / 1e6,
    assessable_genome_mb = A / 1e6,
    coverage_fraction = genome / 1e6 / autosome_mb,
    assessable_fraction = share(A, genome),
    loh_neutral_small_events = neutral[["small"]],
    loh_loss_small_events = loss[["small"]],
    complete_loh_fraction = afrac(one),
    complete_loh_loss_fraction = afrac(one & is_loss),
    complete_loh_neutral_fraction = afrac(one & is_neutral),
    complete_loh_gain_fraction = afrac(one & is_gain),
    relative_loss_fraction = afrac(is_loss),
    no_homologue_fraction = afrac(none),
    neutral_fraction = afrac(is_neutral),
    loh_copy_neutral_share = share(afrac(one & is_neutral), afrac(one)),
    # homozygous deletion size (genome lost entirely; large values suggest a
    # misfit solution rather than biology)
    homdel_mb = sum(len[state == "homdel"], na.rm = TRUE) / 1e6,
    # largest run of consecutive homozygously deleted segments (deleted
    # sequence only; not merged across gaps like the events)
    largest_homdel_mb = {
      hd = !is.na(state) & state == "homdel"
      run = cumsum(c(TRUE, !hd[-1] | !hd[-length(hd)] | x$chr[-1] != x$chr[-nrow(x)]))
      if(any(hd)) max(tapply(len[hd], run[hd], sum)) / 1e6 else 0
    },
    stringsAsFactors = FALSE)
  summary = cbind(summary, extra)

  # copy-neutral LOH events, by extent, under the chosen reference
  cnloh = .battenberg_cnloh(x, a, assessable, len, sample, arms, loci,
                            purity = if(!is.null(cp)) cp$cellularity[1] else NA_real_,
                            doubled = doubled, reference = cnloh_reference,
                            min_event_size = min_event_size,
                            merge_max_gap = merge_max_gap,
                            arm_fraction = arm_fraction, assessable_bp = A)
  summary = cbind(summary, cnloh$summary)

  # genome in each allele-specific state (e.g. 1+1 vs 2+0 vs 2+1)
  allelic_states = NULL
  if(A > 0){
    key = a$allelic_state[assessable]
    mb = tapply(len[assessable], key, sum)
    first = match(names(mb), key)
    allelic_states = data.frame(
      sample = sample, allelic_state = names(mb),
      major_cn = a$major_cn[assessable][first],
      minor_cn = a$minor_cn[assessable][first],
      total_cn = a$total_cn[assessable][first],
      cn_change = a$cn_change[assessable][first],
      homologues = a$homologues[assessable][first],
      balance = a$balance[assessable][first],
      genome_mb = as.numeric(mb) / 1e6,
      assessable_fraction = as.numeric(mb) / A,
      stringsAsFactors = FALSE)
    allelic_states = allelic_states[order(allelic_states$total_cn,
                                          -allelic_states$major_cn), ]
    rownames(allelic_states) = NULL
  }

  # every autosomal segment with its observed BAF / LogR-based total and
  # clonal allele state, e.g. for plotting a sample's profile
  segments = data.frame(sample = sample, chrom = x$chr, start = x$startpos,
                        end = x$endpos, BAF = x$BAF, LogR = x$LogR,
                        ntot = x$ntot, major = maj, minor = min, total = tot,
                        subclonal = !is.na(x$frac1_A) & x$frac1_A < 1,
                        state = state,
                        major_cn = a$major_cn, minor_cn = a$minor_cn,
                        total_cn = a$total_cn, allelic_state = a$allelic_state,
                        cn_change = a$cn_change, homologues = a$homologues,
                        balance = a$balance, assessable = assessable,
                        stringsAsFactors = FALSE)
  list(summary = summary, loh_segments = events, segments = segments,
       allelic_states = allelic_states, complete_loh_events = cnloh$events)
}

# Merged events for one sample: runs of consecutive segments with the same
# (non-NA) class, merged across interruptions of up to merge_max_gap, with
# size, terminal flag, arm and arm-level flag. NULL if there are none.
.battenberg_events = function(x, event_class, sample, major, total,
                              merge_max_gap, arms, arm_fraction,
                              first_pos, last_pos, merge_uncalled = FALSE){
  keep = !is.na(event_class)
  if(!any(keep)){
    return(NULL)
  }
  run_id = cumsum(c(TRUE, x$chr[-1] != x$chr[-nrow(x)] |
                      is.na(event_class[-1]) | is.na(event_class[-nrow(x)]) |
                      event_class[-1] != event_class[-nrow(x)]))
  events = do.call(rbind, lapply(split(which(keep), run_id[keep]), function(ix){
    chrom = x$chr[ix[1]]
    start = min(x$startpos[ix]); end = max(x$endpos[ix])
    data.frame(sample = sample, chrom = chrom, start = start, end = end,
               size = end - start, class = event_class[ix[1]],
               major = paste(unique(major[ix]), collapse = ","),
               total_cn = paste(unique(total[ix]), collapse = ","),
               stringsAsFactors = FALSE)
  }))
  # merge events of the same class separated by short interruptions, of
  # any state (e.g. a homozygous deletion inside a copy-neutral LOH region);
  # interrupting events are kept as events of their own
  events = events[order(as.integer(events$chrom), events$start), ]
  combine = function(a, b) paste(unique(c(strsplit(a, ",")[[1]],
                                          strsplit(b, ",")[[1]])),
                                 collapse = ",")
  merged = list()
  for(i in seq_len(nrow(events))){
    e = events[i, ]
    same = which(vapply(merged, function(m) m$chrom == e$chrom &&
                          m$class == e$class, logical(1)))
    target = if(length(same) > 0) same[length(same)] else NA
    # the interruption: the whole gap, or (merge_uncalled) only the segments
    # in it, so uncalled sequence such as centromeres doesn't separate events
    gap = if(!is.na(target)){
      if(merge_uncalled){
        between = x$chr == e$chrom & x$startpos >= merged[[target]]$end &
          x$endpos <= e$start
        sum(as.numeric(x$endpos[between] - x$startpos[between]))
      }else e$start - merged[[target]]$end
    }else Inf
    if(!is.na(target) && gap <= merge_max_gap){
      merged[[target]]$end = max(merged[[target]]$end, e$end)
      merged[[target]]$major = combine(merged[[target]]$major, e$major)
      merged[[target]]$total_cn = combine(merged[[target]]$total_cn, e$total_cn)
    }else{
      merged[[length(merged) + 1]] = e
    }
  }
  events = do.call(rbind, merged)
  events$size = events$end - events$start
  events$terminal = events$start <= first_pos[events$chrom] |
    events$end >= last_pos[events$chrom]
  rownames(events) = NULL
  # arm(s) overlapped and whether the event covers most of an arm
  events$arm = NA_character_
  events$arm_level = FALSE
  for(i in seq_len(nrow(events))){
    a = arms[arms$chromosome == events$chrom[i] &
               arms$start < events$end[i] & arms$end > events$start[i], ]
    if(nrow(a) == 0) next
    covered = pmin(a$end, events$end[i]) - pmax(a$start, events$start[i])
    events$arm[i] = paste0(events$chrom[i], a$arm, collapse = ",")
    events$arm_level[i] = any(covered / (a$end - a$start) >= arm_fraction)
  }
  events
}

# Allele-specific properties of each segment, robust to subclonality.
# Allele copy numbers must be integers (Battenberg's are); other values are
# treated as unknown rather than rounded. For a subclonal segment (two
# states), a property is taken from the dominant state if it covers at least
# min_clonal_fraction of tumour cells, otherwise only when both states agree
# on it, and is NA when they disagree. The baseline is the modal total copy
# number (by length) over segments whose total is known this way.
.battenberg_allele_properties = function(x, min_clonal_fraction){
  as_int = function(v){
    v = suppressWarnings(as.numeric(v))
    ifelse(!is.na(v) & abs(v - round(v)) < 1e-6, round(v), NA_real_)
  }
  maj1 = as_int(x$nMaj1_A); min1 = as_int(x$nMin1_A)
  maj2 = as_int(x$nMaj2_A); min2 = as_int(x$nMin2_A)
  f1 = suppressWarnings(as.numeric(x$frac1_A))
  f2 = suppressWarnings(as.numeric(x$frac2_A))
  # a second state counts only if present and non-zero
  has2 = !is.na(f2) & f2 > 0 & !is.na(maj2) & !is.na(min2)
  clonal = !has2
  dominant_is_2 = has2 & !is.na(f1) & f2 > f1
  dominant_frac = ifelse(dominant_is_2, f2, f1)
  confident = clonal | (!is.na(dominant_frac) & dominant_frac >= min_clonal_fraction)
  majd = ifelse(dominant_is_2, maj2, maj1); mind = ifelse(dominant_is_2, min2, min1)
  pick = function(v1, v2, vd){
    ifelse(clonal, v1, ifelse(confident, vd, ifelse(!is.na(v1) & !is.na(v2) & v1 == v2, v1, NA)))
  }
  tot1 = maj1 + min1; tot2 = maj2 + min2; totd = majd + mind
  total = pick(tot1, tot2, totd)
  len = as.numeric(x$endpos - x$startpos)
  known = !is.na(total)
  baseline = if(any(known)) as.numeric(names(which.max(tapply(len[known], total[known], sum)))) else NA_real_
  change_of = function(t) ifelse(is.na(t), NA, ifelse(t < baseline, "loss", ifelse(t > baseline, "gain", "neutral")))
  homologues_of = function(t, b) ifelse(is.na(t) | is.na(b), NA, ifelse(t == 0, "none", ifelse(b == 0, "one", "both")))
  balance_of = function(a, b) ifelse(is.na(a) | is.na(b) | b == 0, NA, ifelse(a == b, "balanced", "imbalanced"))
  major_cn = pick(maj1, maj2, majd)
  minor_cn = pick(min1, min2, mind)
  # exact allele copy numbers only when the state itself is determined
  exact = clonal | confident | (!is.na(maj1) & !is.na(maj2) & maj1 == maj2 & min1 == min2)
  major_cn[!exact] = NA; minor_cn[!exact] = NA
  list(baseline = baseline,
       major_cn = major_cn, minor_cn = minor_cn,
       total_cn = total,
       allelic_state = ifelse(is.na(major_cn) | is.na(minor_cn), NA_character_,
                              paste0(major_cn, "+", minor_cn)),
       cn_change = pick(change_of(tot1), change_of(tot2), change_of(totd)),
       homologues = pick(homologues_of(tot1, min1), homologues_of(tot2, min2), homologues_of(totd, mind)),
       balance = pick(balance_of(maj1, min1), balance_of(maj2, min2), balance_of(majd, mind)))
}

# Complete LOH events for one sample, built from assessable, complete LOH
# segments (one homologue, total > 0) with an exact total copy number:
# consecutive segments with the same total are merged (across uncalled
# sequence and up to merge_max_gap of other segments). Each event's loh_type
# (neutral, loss or gain) is relative to the chosen `reference`. Each event is classified by extent,
# relative to the called part of its chromosome/arm: whole_chromosome
# (covers at least arm_fraction of the chromosome's called span), arm (at
# least arm_fraction of an arm's called span), terminal (reaches the first or
# last called position) or interstitial. Returns the events table and the
# per-sample summary columns for the chosen reference.
.battenberg_cnloh = function(x, a, assessable, len, sample, arms, loci,
                             purity, doubled, reference, min_event_size,
                             merge_max_gap, arm_fraction, assessable_bp){
  baseline = a$baseline
  total = a$total_cn
  candidate = assessable & !is.na(a$homologues) & a$homologues == "one" &
    !is.na(total) & total > 0 & !is.na(baseline)
  counts_as = function(t) if(reference == "sample") t == baseline else t == 2
  first_pos = tapply(x$startpos, x$chr, min)
  last_pos = tapply(x$endpos, x$chr, max)
  cls = ifelse(candidate, paste0("loh_", total), NA_character_)
  events = .battenberg_events(x, cls, sample, major = a$major_cn, total = total,
                              merge_max_gap = merge_max_gap, arms = arms,
                              arm_fraction = arm_fraction,
                              first_pos = first_pos, last_pos = last_pos,
                              merge_uncalled = TRUE)
  # the event each candidate segment belongs to
  seg_event = rep(NA_integer_, nrow(x))
  if(!is.null(events)){
    for(i in which(candidate)){
      hit = which(events$chrom == x$chr[i] & events$class == cls[i] &
                    events$start <= x$startpos[i] & events$end >= x$endpos[i])
      if(length(hit) > 0) seg_event[i] = hit[1]
    }
    events$total_cn = as.numeric(sub("^loh_", "", events$class))
    events$class = NULL
    # coverage of the called span of the chromosome and of each arm
    events$extent = NA_character_
    events$arm_level = FALSE
    for(i in seq_len(nrow(events))){
      chrom = events$chrom[i]
      span = c(first_pos[[chrom]], last_pos[[chrom]])
      covered = function(lo, hi) max(0, min(hi, events$end[i]) - max(lo, events$start[i])) / (hi - lo)
      chrom_arms = arms[arms$chromosome == chrom, , drop = FALSE]
      arm_cov = vapply(seq_len(nrow(chrom_arms)), function(k){
        lo = max(chrom_arms$start[k], span[1]); hi = min(chrom_arms$end[k], span[2])
        if(hi <= lo) 0 else covered(lo, hi)
      }, numeric(1))
      events$arm_level[i] = any(arm_cov >= arm_fraction)
      events$extent[i] = if(covered(span[1], span[2]) >= arm_fraction) "whole_chromosome"
        else if(events$arm_level[i]) "arm"
        else if(events$terminal[i]) "terminal" else "interstitial"
    }
    events$allelic_state = vapply(seq_len(nrow(events)), function(i){
      s = a$allelic_state[which(seg_event == i)]
      paste(unique(stats::na.omit(s)), collapse = ",")
    }, character(1))
    events$baseline_cn = baseline
    events$doubled = doubled
    events$neutral_vs_sample = events$total_cn == baseline
    events$neutral_vs_diploid = events$total_cn == 2
    ref_total = if(reference == "sample") baseline else 2
    events$loh_type = ifelse(events$total_cn == ref_total, "neutral",
                             ifelse(events$total_cn < ref_total, "loss", "gain"))
    events$large = events$size >= min_event_size
    events$counted = events$loh_type == "neutral" & events$large
    # share of the event's LOH sequence in subclonal segments
    subclonal = !is.na(x$frac1_A) & x$frac1_A < 1
    events$subclonal_fraction = vapply(seq_len(nrow(events)), function(i){
      ix = which(seg_event == i)
      sum(len[ix] * subclonal[ix]) / sum(len[ix])
    }, numeric(1))
    # length-weighted mean |observed - expected| major-allele BAF at the
    # fitted purity, for a + 0 states
    baf_major = pmax(x$BAF, 1 - x$BAF)
    expected = if(is.na(purity)) rep(NA_real_, nrow(x)) else
      (purity * a$major_cn + (1 - purity)) / (purity * a$major_cn + 2 * (1 - purity))
    events$baf_deviation = vapply(seq_len(nrow(events)), function(i){
      ix = which(seg_event == i & !is.na(expected) & !is.na(baf_major))
      if(length(ix) == 0) NA_real_ else sum(len[ix] * abs(baf_major[ix] - expected[ix])) / sum(len[ix])
    }, numeric(1))
    events$genes = vapply(seq_len(nrow(events)), function(i){
      if(is.null(loci)) return("")
      g = loci$gene[loci$chrom == events$chrom[i] & loci$start < events$end[i] &
                      loci$end > events$start[i]]
      paste(unique(stats::na.omit(g)), collapse = ",")
    }, character(1))
    events = events[, c("sample", "chrom", "start", "end", "size", "arm",
                        "extent", "allelic_state", "total_cn", "baseline_cn",
                        "doubled", "loh_type", "neutral_vs_sample",
                        "neutral_vs_diploid", "large", "counted",
                        "subclonal_fraction", "baf_deviation", "genes")]
  }
  loh_events = events

  # summary for the chosen reference: counts and Mb from events of at least
  # min_event_size; fractions from all segments, of the assessable genome
  counted = if(is.null(events)) events else events[events$counted, , drop = FALSE]
  if(!is.null(counted) && nrow(counted) == 0) counted = NULL
  n_ext = function(e) if(is.null(counted)) 0L else sum(counted$extent == e)
  mb_ext = function(e) if(is.null(counted)) 0 else sum(counted$size[counted$extent == e]) / 1e6
  in_ref = candidate & counts_as(total)
  # segments in whole-chromosome LOH events of any copy number
  whole = rep(FALSE, nrow(x))
  if(!is.null(loh_events)){
    whole[!is.na(seg_event)] = loh_events$extent[seg_event[!is.na(seg_event)]] == "whole_chromosome"
  }
  # segmental LOH events of at least min_event_size, and the copy-neutral ones
  sized = if(is.null(loh_events)) loh_events else
    loh_events[loh_events$size >= min_event_size & loh_events$extent != "whole_chromosome", ]
  n_sized = if(is.null(sized)) 0 else nrow(sized)
  segmental = if(is.null(counted)) counted else counted[counted$extent != "whole_chromosome", ]
  summary = data.frame(
    cnloh_reference = reference,
    cnloh_events = if(is.null(counted)) 0L else nrow(counted),
    cnloh_whole_chromosome_events = n_ext("whole_chromosome"),
    cnloh_arm_events = n_ext("arm"),
    cnloh_terminal_events = n_ext("terminal"),
    cnloh_interstitial_events = n_ext("interstitial"),
    cnloh_whole_chromosome_mb = mb_ext("whole_chromosome"),
    cnloh_arm_mb = mb_ext("arm"),
    cnloh_terminal_mb = mb_ext("terminal"),
    cnloh_interstitial_mb = mb_ext("interstitial"),
    cnloh_largest_segmental_mb = if(is.null(segmental) || nrow(segmental) == 0) 0
      else max(segmental$size) / 1e6,
    cnloh_fraction = if(assessable_bp > 0) sum(len[in_ref]) / assessable_bp else NA_real_,
    cnloh_segmental_fraction = if(assessable_bp > 0)
      sum(len[in_ref & !whole]) / assessable_bp else NA_real_,
    loh_copy_neutral_share_segmental = {
      d = sum(len[candidate & !whole])
      if(d > 0) sum(len[in_ref & !whole]) / d else NA_real_
    },
    loh_copy_neutral_event_share = if(n_sized > 0)
      sum(counts_as(sized$total_cn)) / n_sized else NA_real_,
    stringsAsFactors = FALSE)
  list(events = events, summary = summary)
}
