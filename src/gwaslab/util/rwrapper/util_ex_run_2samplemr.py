import os
from typing import TYPE_CHECKING, Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from gwaslab.g_SumstatsPair import SumstatsPair

from gwaslab.extension import _check_susie_version, _checking_r_version
from gwaslab.info.g_Log import Log
from gwaslab.qc.qc_decorator import with_logging
from gwaslab.util.general.util_ex_result_manager import ResultManager
from gwaslab.util.rwrapper.util_ex_r_runner import RScriptRunner
from gwaslab.util.util_in_convert_h2 import _get_per_snp_r2


@with_logging(
    start_to_msg="run MR using twosampleMR from command line",
    finished_msg="running MR using twosampleMR from command line",
    start_cols=["SNPID", "CHR", "POS", "EA", "NEA"],
    start_function=".run_two_sample_mr()",
    check_tools=["r", "r:TwoSampleMR"]
)
def _run_two_sample_mr(
    sumstatspair_object: "SumstatsPair",
    r: str,
    out: str = "./",
    clump: bool = False,
    f_check: float = 10,
    exposure1: str = "Trait1",
    outcome2: str = "Trait2",
    n1: int | None = None,
    n2: int | None = None,
    binary1: bool = False,
    cck1: tuple[int, int, float] | None = None,
    cck2: tuple[int, int, float] | None = None,
    ncase1: int | None = None,
    ncontrol1: int | None = None,
    prevalence1: float | None = None,
    binary2: bool = False,
    ncase2: int | None = None,
    ncontrol2: int | None = None,
    prevalence2: float | None = None,
    methods: list[str] | None = None,
    log: Log = Log()
) -> None:

    if methods is None:
        methods = ["mr_ivw","mr_simple_mode","mr_weighted_median","mr_egger_regression","mr_ivw_mre", "mr_weighted_mode"]
        methods_string = '"{}"'.format('","'.join(methods))

    if cck1 is not None:
        log.write(f" - ncase1, ncontrol1, prevalence1:{cck1}")
        binary1 = True
        ncase1 = cck1[0]
        ncontrol1 = cck1[1]
        prevalence1 =  cck1[2]
        n1 = ncase1 + ncontrol1
    if cck2 is not None:
        log.write(f" - ncase2, ncontrol2, prevalence2:{cck2}")
        binary2 = True
        ncase2 = cck2[0]
        ncontrol2 = cck2[1]
        prevalence2 =  cck2[2]
        n2 = ncase2 + ncontrol2

    if clump==True:
        sumstatspair = sumstatspair_object.clumps["clumps"]
    else:
        sumstatspair = sumstatspair_object.data

    if n1 is not None:
        sumstatspair["N_1"] = n1
    if n2 is not None:
        sumstatspair["N_2"] = n2

    sumstatspair = _filter_by_f(sumstatspair, f_check, n1, binary1, ncase1, prevalence1)

    log = _checking_r_version(r, log)
    #log = _check_susie_version(r,log)

    cols_for_trait1, cols_for_trait2 = _sort_columns_to_load(sumstatspair)
    cols_for_trait1_script = _cols_list_to_r_script(cols_for_trait1)
    cols_for_trait2_script = _cols_list_to_r_script(cols_for_trait2)

    # Initialize runners and managers
    runner = RScriptRunner(
        r=r,
        log=log,
        timeout=None,
        temp_dir=out,
        cleanup=True
    )
    result_manager = ResultManager(log=log)

    # Prepare output paths
    memory_id = id(sumstatspair)
    prefix = f"{exposure1}_{outcome2}_{memory_id}"
    prefix = "{}{}".format(out.rstrip("/") + "/", prefix)
    temp_sumstats_path = "{out}twosample_mr_{exposure}_{outcome}_{memory_id}.csv.gz".format(
        out=out.rstrip("/") + "/",
        exposure=exposure1,
        outcome=outcome2,
        memory_id=memory_id
    )

    if len(sumstatspair) > 0:
        sumstatspair.to_csv(temp_sumstats_path, index=None)
    else:
        return 0
    ###
    calculate_r_script = ""

    if binary1==True:
        calculate_r_script+= _make_script_for_calculating_r("exposure", ncase1, ncontrol1, prevalence1)
    else:
        calculate_r_script+= _make_script_for_calculating_r_quant("exposure")

    if binary2==True:
        calculate_r_script+= _make_script_for_calculating_r("outcome", ncase2, ncontrol2, prevalence2)
    else:
        calculate_r_script+= _make_script_for_calculating_r_quant("outcome")

    # create scripts
    directionality_test_script=f"""
    results_directionality <- directionality_test(harmonized_data)
    write.csv(results_directionality, "{prefix}.directionality", row.names = FALSE)
    """

    # Two Sample MR
    # Tests
    ## Pleiotropy
    ## Heterogeneity
    ## Reverse causality
    rscript="""
    library(TwoSampleMR)
    sumstats <- read.csv("{temp_sumstats_path}")

    {pheno1}
    {pheno2}

    exp_raw <-sumstats[,c({cols_for_trait1})]

    out_raw <-sumstats[,c({cols_for_trait2})] 

    exp_dat <- format_data( exp_raw,
                  type = "exposure",
                  snp_col = "SNPID",
                  beta_col = "BETA_1",
                  se_col = "SE_1",
                  effect_allele_col = "EA",
                  other_allele_col = "NEA",
                  eaf_col = "EAF_1",
                  pval_col = "P_1",
                  phenotype_col = "PHENO_1",
                  samplesize_col= "N_1"
                 )
    
    out_dat <- format_data( out_raw,
                  type = "outcome",
                  snp_col = "SNPID",
                  beta_col = "BETA_2",
                  se_col = "SE_2",
                  effect_allele_col = "EA",
                  other_allele_col = "NEA",
                  eaf_col = "EAF_2",
                  pval_col = "P_2",
                  phenotype_col = "PHENO_2",
                  samplesize_col= "N_2"
                  )
    
    harmonized_data <- harmonise_data(exp_dat,out_dat,action=1)
    
    results_mr <- mr(harmonized_data, method_list = c({methods_string}))
    write.csv(results_mr, "{prefix}.mr", row.names = FALSE)

    results_heterogeneity <- mr_heterogeneity(harmonized_data)
    write.csv(results_heterogeneity, "{prefix}.heterogeneity", row.names = FALSE)

    results_pleiotropy <- mr_pleiotropy_test(harmonized_data)
    write.csv(results_pleiotropy, "{prefix}.pleiotropy", row.names = FALSE)
    
    results_single <- mr_singlesnp(harmonized_data)
    write.csv(results_single, "{prefix}.singlesnp", row.names = FALSE)
    
    results_loo <- mr_leaveoneout(harmonized_data)
    write.csv(results_loo, "{prefix}.leaveoneout", row.names = FALSE)
    
    {calculate_r}
    {directionality_test}

    """.format(
        temp_sumstats_path = temp_sumstats_path,
        pheno1 = f'sumstats$PHENO_1 <- "{exposure1}"' if exposure1 is not None else "",
        pheno2 = f'sumstats$PHENO_2 <- "{outcome2}"' if outcome2 is not None else "",
        cols_for_trait1 = cols_for_trait1_script,
        cols_for_trait2 = cols_for_trait2_script,
        methods_string=methods_string,
        prefix=prefix,
        calculate_r = calculate_r_script,
        directionality_test = directionality_test_script
    )

    # Prepare expected output files
    expected_outputs = [
        f"{prefix}.mr",
        f"{prefix}.pleiotropy",
        f"{prefix}.heterogeneity",
        f"{prefix}.singlesnp",
        f"{prefix}.leaveoneout",
        f"{prefix}.directionality"
    ]

    # Execute R script
    log.write(" Running TwoSampleMR from command line...")
    result = runner.execute(
        script_content=rscript,
        expected_outputs=expected_outputs,
        temp_prefix=f"twosamplemr_{exposure1}_{outcome2}",
        temp_suffix=".R",
        timeout=None,
        verbose=True,
        working_dir=out.rstrip("/") + "/"
    )

    # Trace result
    result_manager.trace(
        result=result,
        identifier=f"{exposure1}_{outcome2}",
        parameters={
            "exposure": exposure1,
            "outcome": outcome2,
            "methods": methods,
            "clump": clump,
            "f_check": f_check
        }
    )

    # Process results
    if result.success:
        # Read output files from result
        for suffix in ["mr", "pleiotropy", "heterogeneity", "singlesnp", "leaveoneout", "directionality"]:
            expected_file = f"{prefix}.{suffix}"
            if expected_file in result.output_files:
                try:
                    sumstatspair_object.mr[suffix] = pd.read_csv(result.output_files[expected_file])
                except Exception as e:
                    log.warning(f"Error reading {suffix} output: {e!s}")

        # Store R log
        sumstatspair_object.mr["r_log"] = result.output
    else:
        log.write(" Error during TwoSampleMR execution!")
        if result.errors:
            for error in result.errors:
                log.warning(f"  - {error}")
        log.write(" R script content:")
        log.write(rscript)
        if result.output:
            log.write(" R output:")
            log.write(result.output)



def _sort_columns_to_load(sumstatspair: pd.DataFrame) -> tuple[list[str], list[str]]:
    cols_for_trait1=["SNPID","CHR","POS","EA","NEA","PHENO_1","N_1"]
    cols_for_trait2=["SNPID","CHR","POS","EA","NEA","PHENO_2","N_2"]
    for i in ["EAF","BETA","SE","P"]:
        if i+"_1" not in cols_for_trait1 and i+"_1" in sumstatspair.columns:
            cols_for_trait1.append(i+"_1")
        if i+"_2" not in cols_for_trait2 and i+"_2" in sumstatspair.columns:
            cols_for_trait2.append(i+"_2")
    return cols_for_trait1, cols_for_trait2




def _cols_list_to_r_script(cols_for_trait1: list[str]) -> str:
    script = '"{}"'.format('","'.join(cols_for_trait1))
    return script

def _make_script_for_calculating_r(
    exposure_or_outcome: str,
    ncase: int,
    ncontrol: int,
    prevalence: float
) -> str:

        script = f"""
        harmonized_data$"r.{exposure_or_outcome}" <- get_r_from_lor(  harmonized_data$"beta.{exposure_or_outcome}",
                                                        harmonized_data$"eaf.{exposure_or_outcome}",
                                                        {ncase},
                                                        {ncontrol},
                                                        {prevalence},
                                                        model = "logit",
                                                        correction = FALSE
                                                        )
        """
        return script


def _make_script_for_calculating_r_quant(exposure_or_outcome: str) -> str:
        script = f"""
        harmonized_data$"r.{exposure_or_outcome}" <- get_r_from_bsen(  harmonized_data$"beta.{exposure_or_outcome}",
                                                        harmonized_data$"se.{exposure_or_outcome}",
                                                        harmonized_data$"samplesize.{exposure_or_outcome}"
                                                        )
        """
        return script


def _filter_by_f(
    sumstatspair: pd.DataFrame,
    f_check: float,
    n1: int | None,
    binary1: bool | None = None,
    ncase1: int | None = None,
    ncontrol1: int | None = None,
    prevalence1: float | None = None,
    log: Log = Log()
) -> pd.DataFrame:

    if binary1==True:
        sumstatspair = _get_per_snp_r2(sumstatspair,
            beta="BETA_1",
            af="EAF_1",
            n = "N_1",
            mode="b",
            se="SE_1",
            vary=1,
            ncase=ncase1,
            ncontrol=ncontrol1,
            prevalence=prevalence1,
            k=1)
    else:
        sumstatspair = _get_per_snp_r2(sumstatspair,
            beta="BETA_1",
            af="EAF_1",
            n = "N_1",
            mode="q",
            se="SE_1",
            vary=1,
            k=1)

    log.write("Filtered out {} variants with F < {}".format(len(sumstatspair) - sum(sumstatspair["F"]>f_check),f_check))
    sumstatspair = sumstatspair.loc[sumstatspair["F"]>f_check,:]

    return sumstatspair
