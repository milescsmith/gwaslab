import copy
from typing import TYPE_CHECKING, Any, List, Optional, Union

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.info.g_Log import Log

import gc

from gwaslab.bd.bd_common_data import (
    get_chr_list,
    get_chr_to_number,
    get_format_dict,
    get_formats_list,
    get_high_ld,
    get_number_to_chr,
)
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap3
from gwaslab.hm.hm_casting import _align_with_mold, _merge_mold_with_sumstats_by_chrpos
from gwaslab.hm.hm_harmonize_sumstats import (
    _check_ref,
    _parallele_infer_af_with_maf,
    _parallelize_assign_rsid,
    _parallelize_check_af,
    _parallelize_infer_af,
    _parallelize_infer_strand,
    _parallelize_rsid_to_chrpos,
)
from gwaslab.hm.hm_liftover_v2 import _liftover_variant
from gwaslab.info.g_Log import Log
from gwaslab.info.g_meta import _append_meta_record, _init_meta, _update_meta
from gwaslab.info.g_Sumstats_summary import lookupstatus, summarize
from gwaslab.info.g_version import _show_version, gwaslab_info
from gwaslab.io.io_load_ld import _to_finemapping_using_ld
from gwaslab.io.io_preformat_input_polars import preformatp
from gwaslab.io.io_read_pipcs import _read_pipcs
from gwaslab.io.io_to_formats import _to_format
from gwaslab.qc.qc_build import _set_build
from gwaslab.qc.qc_fix_sumstats import (
    _fix_allele,
    _fix_chr,
    _fix_ID,
    _fix_pos,
    _flip_allele_stats,
    _flip_SNPID,
    _parallelize_normalize_allele,
    _process_build,
    _remove_dup,
    _sort_column,
    _sort_coordinate,
    _strip_SNPID,
)
from gwaslab.qc.qc_fix_sumstats_polars import _fix_chrp, _fix_posp
from gwaslab.qc.qc_sanity_check import _check_data_consistency, _sanity_check_stats
from gwaslab.util.rwrapper.util_ex_run_susie import _get_cs_lead, _run_susie_rss
from gwaslab.util.util_abf_finemapping import _abf_finemapping, _make_cs
from gwaslab.util.util_ex_calculate_ldmatrix import _to_finemapping
from gwaslab.util.util_ex_calculate_prs import _calculate_prs
from gwaslab.util.util_ex_ldproxyfinder import _extract_ld_proxy
from gwaslab.util.util_ex_ldsc import (
    _estimate_h2_by_ldsc,
    _estimate_h2_cts_by_ldsc,
    _estimate_partitioned_h2_by_ldsc,
    _estimate_rg_by_ldsc,
)
from gwaslab.util.util_ex_run_clumping import _clump
from gwaslab.util.util_ex_run_prscs import _run_prscs
from gwaslab.util.util_in_calculate_gc import _lambda_GC
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2
from gwaslab.util.util_in_estimate_ess import _get_ess
from gwaslab.util.util_in_fill_data import _fill_data
from gwaslab.util.util_in_filter_value import (
    _exclude_hla,
    _filter_in,
    _filter_indel,
    _filter_out,
    _filter_palindromic,
    _filter_region,
    _filter_region_in,
    _filter_region_out,
    _filter_snp,
    _filter_values,
    _get_flanking,
    _get_flanking_by_chrpos,
    _get_flanking_by_id,
    _infer_build,
    _sampling,
    _search_variants,
)
from gwaslab.util.util_in_get_density import _get_signal_density2
from gwaslab.util.util_in_get_sig import _anno_gene, _check_cis, _check_novel_set, _get_novel, _get_sig
from gwaslab.viz.viz_plot_compare_af import plotdaf
from gwaslab.viz.viz_plot_credible_sets import _plot_cs
from gwaslab.viz.viz_plot_mqqplot import _mqqplot
from gwaslab.viz.viz_plot_phe_heatmap import _gwheatmap
from gwaslab.viz.viz_plot_trumpetplot import _plot_trumpet


#20220309
class Sumstatsp:
    def __init__(
        self,
        sumstats: Union[str, pd.DataFrame, Any],
        fmt: str | None = None,
        tab_fmt: str = "tsv",
        snpid: str | None = None,
        rsid: str | None = None,
        chrom: str | None = None,
        pos: str | None = None,
        ea: str | None = None,
        nea: str | None = None,
        ref: str | None = None,
        alt: str | None = None,
        eaf: str | None = None,
        neaf: str | None = None,
        maf: str | None = None,
        n: str | None = None,
        beta: str | None = None,
        se: str | None = None,
        chisq: str | None = None,
        z: str | None = None,
        f: str | None = None,
        t: str | None = None,
        p: str | None = None,
        q: str | None = None,
        mlog10p: str | None = None,
        test: str | None = None,
        info: str | None = None,
        OR: str | None = None,
        OR_95L: str | None = None,
        OR_95U: str | None = None,
        beta_95L: str | None = None,
        beta_95U: str | None = None,
        HR: str | None = None,
        HR_95L: str | None = None,
        HR_95U: str | None = None,
        ncase: str | None = None,
        ncontrol: str | None = None,
        neff: str | None = None,
        i2: str | None = None,
        phet: str | None = None,
        dof: str | None = None,
        snpr2: str | None = None,
        status: str | None = None,
        other: list[str] = [],
        chrom_pat: str | None = None,
        snpid_pat: str | None = None,
        usekeys: list[str] | None = None,
        direction: str | None = None,
        verbose: bool = True,
        study: str = "Study_1",
        trait: str = "Trait_1",
        build: str = "99",
        species: str = "homo sapiens",
        build_infer: bool = False,
        **readargs: Any
    ) -> None:

        # basic attributes
        self.data = pd.DataFrame()
        self.log = Log()
        self.ldsc_h2 = None
        self.ldsc_h2_results = None
        self.ldsc_rg = pd.DataFrame()
        self.ldsc_h2_cts = None
        self.ldsc_partitioned_h2_summary = None
        self.ldsc_partitioned_h2_results = None
        # meta information
        self.meta = _init_meta()
        self.build = build
        self.meta["gwaslab"]["study_name"] =  study
        self.meta["gwaslab"]["species"] = species

        # initialize attributes for clumping and finmapping
        #self.to_finemapping_file_path = ""
        #self.to_finemapping_file  = pd.DataFrame()
        #self.plink_log = ""

        # path / file / plink_log
        self.finemapping = dict()

        # clumps / clumps_raw / plink_log
        self.clumps = dict()

        self.pipcs = pd.DataFrame()

        # print gwaslab version information
        _show_version(self.log, verbose=verbose)

        #preformat the data
        self.data  = preformatp(
          sumstats=sumstats,
          fmt=fmt,
          tab_fmt = tab_fmt,
          snpid=snpid,
          rsid=rsid,
          chrom=chrom,
          pos=pos,
          ea=ea,
          nea=nea,
          ref=ref,
          alt=alt,
          eaf=eaf,
          neaf=neaf,
          maf=maf,
          n=n,
          beta=beta,
          se=se,
          chisq=chisq,
          z=z,
          f=f,
          t=t,
          p=p,
          q=q,
          mlog10p=mlog10p,
          test=test,
          info=info,
          OR=OR,
          OR_95L=OR_95L,
          OR_95U=OR_95U,
          beta_95L=beta_95L,
          beta_95U=beta_95U,
          HR=HR,
          HR_95L=HR_95L,
          HR_95U=HR_95U,
          i2=i2,
          phet=phet,
          dof=dof,
          snpr2=snpr2,
          ncase=ncase,
          ncontrol=ncontrol,
          neff=neff,
          direction=direction,
          study=study,
          build=build,
          trait=trait,
          status=status,
          other=other,
          usekeys=usekeys,
          chrom_pat=chrom_pat,
          snpid_pat=snpid_pat,
          verbose=verbose,
          readargs=readargs,
          log=self.log)

        gc.collect()

    def fix_chr(self, **kwargs: Any) -> "Sumstatsp":
        """
        Standardize chromosome notation and handle special chromosome cases (X, Y, MT).
        
        This method normalizes chromosome labels to a consistent format, extracts chromosome
        numbers from various formats (e.g., "chr1", "1", "chrX"), maps special chromosomes
        (X, Y, mitochondrial) to standardized numeric identifiers, and optionally removes
        invalid chromosome values.
        
        Parameters
        ----------
        chrom : str, default "CHR"
            Column name for chromosome.
        status : str, default "STATUS"
            Column name for status.
        add_prefix : str, optional, default=""
            Prefix to prepend to chromosome labels (e.g., "chr").
        remove : bool, default False
            If True, remove records with invalid or unrecognized chromosome labels.
        verbose : bool, default True
            If True, print progress or diagnostic messages.
        
        Returns
        -------
        Sumstatsp
            Returns self for method chaining.
        """
        from gwaslab.io.io_process_kwargs import remove_overlapping_kwargs
        kwargs = remove_overlapping_kwargs(kwargs, {"log"})
        self.data = _fix_chrp(self, log=self.log, **kwargs)
        return self

    def fix_pos(self, **kwargs: Any) -> "Sumstatsp":
        """
        Standardize and validate genomic base-pair positions.
        
        This method checks that reported genomic positions fall within valid chromosomal bounds
        and optionally removes invalid entries. It handles string-formatted positions with
        thousands separators, converts positions to Int64 type, and filters out positions
        outside the specified range.
        
        Parameters
        ----------
        pos : str, default "POS"
            Column name for position.
        status : str, default "STATUS"
            Column name for status.
        remove : bool, default False
            If True, remove records with invalid or out-of-range positions.
        verbose : bool, default True
            If True, print progress or diagnostic messages.
        lower_limit : int, optional
            Minimum acceptable genomic position. Default is 0.
        upper_limit : int, optional
            Maximum acceptable genomic position.
        limit : int, default 250000000
            Default upper limit applied when `upper_limit` is not provided.
        
        Returns
        -------
        Sumstatsp
            Returns self for method chaining.
        """
        from gwaslab.io.io_process_kwargs import remove_overlapping_kwargs
        kwargs = remove_overlapping_kwargs(kwargs, {"log"})
        self.data = _fix_posp(self, log=self.log, **kwargs)
        return self
