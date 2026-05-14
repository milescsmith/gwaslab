import os
import subprocess
from typing import TYPE_CHECKING, List, Optional, Tuple

import numpy as np
import pandas as pd

from gwaslab.extension import _checking_plink_version
from gwaslab.info.g_Log import Log

if TYPE_CHECKING:
    pass

def _process_plink_input_files(chrlist: list[int],
                               bfile: str | None = None,
                               pfile: str | None = None,
                               vcf: str | None = None,
                               bgen: str | None = None,
                               sample: str | None = None,
                               threads: int = 1,
                               plink_log: str = "",
                               log: Log = Log(),
                               overwrite: bool = False,
                               bgen_mode: str = "ref-first",
                               convert: str = "bfile",
                               memory: int | None = None,
                               load_bim: bool = False,
                               plink: str = "plink",
                               plink2: str = "plink2") -> tuple[str, str, list[pd.DataFrame], str]:
    """
    Process input files (bfile,pfile,vcf,bgen) to either PLINK1 bed/bim/fam or PLINK2 pgen/psam/pvar. 
    
    Parameters:
    -----------
    load_bim : bool, default=False
        If True, load BIM/PVAR variant information into ref_bims list.
        If False (default), ref_bims will be an empty list.
        **Important**: Pass load_bim=True to populate ref_bims when using VCF or BGEN inputs.
    
    Returns:
    --------
    ref_file_prefix : prefix for either bfile or pfile.
    plink_log : if plink was used, return the log. Otherwise, return an empty string.
    ref_bims : if load_bim is True, return bim files as a list of pd.DataFrame. Otherwise, empty list.
    filetype : either bfile or pfile.
    """

    # Step 1: Initialize list to store BIM/PVAR dataframes if load_bim is True
    ref_bims = []

    # Step 2: Determine the input file type and check if it uses wildcard pattern (@)
    # Returns: filetype (bfile/pfile/vcf/bgen), ref_file_prefix, and is_wild_card flag
    # Note: bfile and pfile require no conversion, while vcf and bgen need to be converted
    # File prefix can be single file or wildcard pattern with @ for chromosome-specific files
    filetype, ref_file_prefix, is_wild_card = _check_file_type(bfile,  pfile,  vcf,  bgen, sample)

    # Step 3: Process bfile format (PLINK1 bed/bim/fam)
    # No conversion needed, just validate files exist and optionally load BIM data
    if filetype == "bfile":
        ref_file_prefix, ref_bims = _process_bfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix,
                                                   ref_bims=ref_bims,
                                                   is_wild_card=is_wild_card,
                                                   log=log,
                                                   load_bim=load_bim)

    # Step 4: Process pfile format (PLINK2 pgen/pvar/psam)
    # No conversion needed, just validate files exist and optionally load PVAR data
    elif filetype == "pfile":
        ref_file_prefix, ref_bims = _process_pfile(chrlist=chrlist,
                                                   ref_file_prefix=ref_file_prefix,
                                                   ref_bims=ref_bims,
                                                   is_wild_card=is_wild_card,
                                                   log=log,
                                                   load_bim=load_bim)

    # Step 5: Process VCF format - convert to bfile or pfile using PLINK2
    # VCF files need conversion, so filetype is updated to the target format (convert parameter)
    elif filetype == "vcf":
        ref_file_prefix, plink_log, ref_bims = _process_vcf(ref_file_prefix=ref_file_prefix,
                                                            chrlist=chrlist,
                                                            ref_bims=ref_bims,
                                                            is_wild_card=is_wild_card,
                                                            log=log,
                                                            plink_log=plink_log,
                                                            threads=threads,
                                                            convert=convert,
                                                            memory=memory,
                                                            overwrite=overwrite,
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        # Update filetype to the converted format (bfile or pfile)
        filetype = convert

    # Step 6: Process BGEN format - convert to bfile or pfile using PLINK2
    # BGEN files need conversion, so filetype is updated to the target format (convert parameter)
    elif filetype == "bgen":
        ref_file_prefix, plink_log, ref_bims = _process_bgen(ref_file_prefix=ref_file_prefix,
                                                            chrlist=chrlist,
                                                            bgen_mode=bgen_mode,
                                                            sample=sample,
                                                            ref_bims=ref_bims,
                                                            is_wild_card=is_wild_card,
                                                            log=log,
                                                            plink_log=plink_log,
                                                            threads=threads,
                                                            convert=convert,
                                                            memory=memory,
                                                            overwrite=overwrite,
                                                            load_bim=load_bim,
                                                            plink=plink,
                                                            plink2=plink2)
        # Update filetype to the converted format (bfile or pfile)
        filetype = convert

    # Step 7: Return processed file prefix, PLINK log (if conversion occurred),
    #         list of BIM/PVAR dataframes (if load_bim=True), and final filetype
    return ref_file_prefix, plink_log, ref_bims, filetype

def _load_single_bim_to_ref_bims(bpfile_prefix: str, ref_bims: list[pd.DataFrame], log: Log) -> list[pd.DataFrame]:
    bim_path =bpfile_prefix+".bim"
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0,1,3,4,5],
                             header=None,
                             dtype={1:"string",0:"category", 3:"int", 4:"string", 5:"string"}).rename(columns={1:"SNPID",0:"CHR_bim",3:"POS_bim",4:"EA_bim",5:"NEA_bim"})
    log.write(f"   -Variants in ref file: {len(single_bim)}")
    ref_bims.append(single_bim)
    return ref_bims

def _load_single_pvar_to_ref_bims(bpfile_prefix: str, ref_bims: list[pd.DataFrame], log: Log) -> list[pd.DataFrame]:
    if os.path.exists(bpfile_prefix+".pvar"):
        bim_path =bpfile_prefix+".pvar"
    elif os.path.exists(bpfile_prefix+".pvar.zst"):
        bim_path =bpfile_prefix+".pvar.zst"
    single_bim = pd.read_csv(bim_path,
                             sep=r"\s+",
                             usecols=[0,1,2,3,4],
                             header=None,
                             comment="#",
                             dtype={2:"string",0:"category", 1:"int", 3:"string", 4:"string"}).rename(columns={2:"SNPID",0:"CHR_bim",1:"POS_bim",3:"EA_bim",4:"NEA_bim"})
    log.write(f"   -Variants in ref file: {len(single_bim)}")
    ref_bims.append(single_bim)
    return ref_bims

def _check_file_type(bfile: str | None = None,
                     pfile: str | None = None,
                     vcf: str | None = None,
                     bgen: str | None = None,
                     sample: str | None = None) -> tuple[str, str, bool]:

    is_wild_card = False
    if bfile is None:
        if pfile is None:
            if vcf is None:
                if (bgen is None) or (sample is None):
                    raise ValueError("You need to provide one from bfile, pfile, bgen, vcf; for bgen, sample file is required.")
                else:
                    if "@" in bgen:
                        is_wild_card = True
                    return "bgen",  bgen, is_wild_card
            else:
                if "@" in vcf:
                    is_wild_card = True
                return "vcf",  vcf, is_wild_card
        else:
            if "@" in pfile:
                is_wild_card = True
            return "pfile",  pfile, is_wild_card
    else:
        if "@" in bfile:
            is_wild_card = True
        return "bfile",  bfile, is_wild_card

def _process_bfile(chrlist: list[int], ref_file_prefix: str, ref_bims: list[pd.DataFrame], is_wild_card: bool, log: Log, load_bim: bool = False) -> tuple[str, list[pd.DataFrame]]:
    if is_wild_card==False:
        is_bim_exist = os.path.exists(ref_file_prefix+".bim")
        is_bed_exist = os.path.exists(ref_file_prefix+".bed")
        is_fam_exist = os.path.exists(ref_file_prefix+".fam")
        if not (is_bim_exist and is_bed_exist  and is_fam_exist):
            raise ValueError(f"PLINK bfiles are missing : {ref_file_prefix} ")
        if load_bim==True:
            ref_bims = _load_single_bim_to_ref_bims(ref_file_prefix, ref_bims, log)
        log.write(f" -Single PLINK bfile as LD reference panel: {ref_file_prefix}")

    else:
        for i in chrlist:
            single_chr_ref_file_prefix = ref_file_prefix.replace("@",str(i))
            is_bim_exist = os.path.exists(single_chr_ref_file_prefix+".bim")
            is_bed_exist = os.path.exists(single_chr_ref_file_prefix+".bed")
            is_fam_exist = os.path.exists(single_chr_ref_file_prefix+".fam")
            if not (is_bim_exist and is_bed_exist  and  is_fam_exist):
                raise ValueError(f"PLINK bfiles for CHR {i} are missing...")
            if load_bim==True:
                ref_bims = _load_single_bim_to_ref_bims(single_chr_ref_file_prefix, ref_bims, log)
        log.write(f" -Split PLINK bfiles for each CHR as LD reference panel: {ref_file_prefix}")

    return ref_file_prefix, ref_bims

def _process_pfile(chrlist: list[int], ref_file_prefix: str, ref_bims: list[pd.DataFrame], is_wild_card: bool, log: Log, load_bim: bool = False) -> tuple[str, list[pd.DataFrame]]:
    if is_wild_card==False:
        is_bim_exist = os.path.exists(ref_file_prefix+".pvar") or os.path.exists(ref_file_prefix+".pvar.zst")
        is_bed_exist = os.path.exists(ref_file_prefix+".pgen")
        is_fam_exist = os.path.exists(ref_file_prefix+".psam")
        if not (is_bim_exist and is_bed_exist  and is_fam_exist):
            raise ValueError(f"PLINK pfiles are missing : {ref_file_prefix} ")
        if load_bim==True:
            ref_bims = _load_single_pvar_to_ref_bims(ref_file_prefix, ref_bims, log)
        log.write(f" -Single PLINK bfile as LD reference panel: {ref_file_prefix}")
    else:
        for i in chrlist:
            single_chr_ref_file_prefix = ref_file_prefix.replace("@",str(i))
            is_bim_exist = os.path.exists(single_chr_ref_file_prefix+".pvar") or os.path.exists(single_chr_ref_file_prefix+".pvar.zst")
            is_bed_exist = os.path.exists(single_chr_ref_file_prefix+".pgen")
            is_fam_exist = os.path.exists(single_chr_ref_file_prefix+".psam")
            if not (is_bim_exist and is_bed_exist  and  is_fam_exist):
                raise ValueError(f"PLINK pfiles for CHR {i} are missing...")
            if load_bim==True:
                ref_bims = _load_single_pvar_to_ref_bims(ref_file_prefix, ref_bims, log)
        log.write(f" -Split PLINK pfiles for each CHR as LD reference panel: {ref_file_prefix}")

    return ref_file_prefix, ref_bims

def _process_vcf(ref_file_prefix: str,
                 chrlist: list[int],
                 ref_bims: list[pd.DataFrame],
                 is_wild_card: bool,
                 log: Log,
                 plink_log: str,
                 threads: int = 1,
                 convert: str = "bfile",
                 memory: int | None = None,
                 overwrite: bool = False,
                 load_bim: bool = False,
                 plink: str = "plink",
                 plink2: str = "plink2") -> tuple[str, str, list[pd.DataFrame]]:
    log.write(f" -Processing VCF : {ref_file_prefix}...")

    #check plink version
    log = _checking_plink_version(plink2=plink2,log=log)

    # file path prefix to return
    if is_wild_card==True:
        ref_file_prefix_converted = ref_file_prefix.replace(".vcf.gz","")
    else:
        ref_file_prefix_converted = ref_file_prefix.replace(".vcf.gz","") + ".@"

    for i in chrlist:
    # for each chr
        log.write(f"  -Processing VCF for CHR {i}...")

        # if multiple bfiles
        if is_wild_card==True:
            vcf_to_load = ref_file_prefix.replace("@",str(i))
             # output bpfile prefix for each chr : chr to chr
            bpfile_prefix = vcf_to_load.replace(".vcf.gz","")
        else:
            vcf_to_load = ref_file_prefix
             # output bpfile prefix for each chr : all to chr
            bpfile_prefix = vcf_to_load.replace(".vcf.gz","") + f".{i}"


        #check if the file exists
        if convert=="bfile":
            is_file_exist = os.path.exists(bpfile_prefix+".bed")
            make_flag = "--make-bed"
        else:
            is_file_exist = os.path.exists(bpfile_prefix+".pgen")
            make_flag = "--make-pgen vzs"

        # figure memory
        if memory is not None:
            memory_flag = f"--memory {memory}"
        else:
            memory_flag = ""

        #if not existing or overwrite is True
        conversion_success = True
        if (not is_file_exist) or overwrite:
            script_vcf_to_bfile = f"""
            {plink2} \
                --vcf {vcf_to_load} \
                --chr {i} \
                {make_flag} \
                --rm-dup force-first \
                --threads {threads}{memory_flag}\
                --out {bpfile_prefix}
            """

            # execute conversion
            try:
                if convert=="bfile":
                    log.write(f"  -Converting VCF to bfile: {bpfile_prefix}.bim/bed/fam...")
                else:
                    log.write(f"  -Converting VCF to pfile: {bpfile_prefix}.pgen/pvar/psam...")
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
                # Verify conversion succeeded by checking if output file exists
                if convert=="bfile":
                    conversion_success = os.path.exists(bpfile_prefix+".bed")
                else:
                    conversion_success = os.path.exists(bpfile_prefix+".pgen")
                # Read PLINK log file if it exists
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file) as f:
                        plink_log += f.read() + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)
                conversion_success = False
                # Try to read PLINK log file even on error
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file) as f:
                        plink_log += f.read() + "\n"
        else:
            log.write(f"  -Plink {convert} for CHR {i} exists: {bpfile_prefix}. Skipping...")

        # Load BIM/PVAR data if requested and conversion was successful
        # Note: load_bim must be True to populate ref_bims, otherwise it remains empty
        if load_bim == True:
            # Only load if file exists (either from successful conversion or pre-existing)
            if convert == "bfile":
                if os.path.exists(bpfile_prefix+".bim"):
                    ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write(f"  -Warning: BIM file not found for CHR {i}: {bpfile_prefix}.bim. Skipping load.")
            elif os.path.exists(bpfile_prefix+".pvar") or os.path.exists(bpfile_prefix+".pvar.zst"):
                ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
            else:
                log.write(f"  -Warning: PVAR file not found for CHR {i}: {bpfile_prefix}.pvar. Skipping load.")
    return ref_file_prefix_converted, plink_log, ref_bims

def _process_bgen(ref_file_prefix: str,
                  chrlist: list[int],
                  ref_bims: list[pd.DataFrame],
                  is_wild_card: bool,
                  log: Log = Log(),
                  plink_log: str = "",
                  sample: str | None = None,
                  bgen_mode: str = "ref-first",
                  threads: int = 1,
                  convert: str = "bfile",
                  memory: int | None = None,
                  overwrite: bool = False,
                  load_bim: bool = False,
                  plink: str = "plink",
                  plink2: str = "plink2") -> tuple[str, str, list[pd.DataFrame]]:
    log.write(f" -Processing BGEN files : {ref_file_prefix}...")

    #check plink version
    log = _checking_plink_version(log=log,plink2=plink2)

    # file path prefix to return
    if is_wild_card==True:
        ref_file_prefix_converted = ref_file_prefix.replace(".bgen","")
    else:
        ref_file_prefix_converted = ref_file_prefix.replace(".bgen","") + ".@"

    for i in chrlist:
    # for each chr
        log.write(f"  -Processing BGEN for CHR {i}...")

        # if multiple bfiles
        if is_wild_card==True:
            bgen_to_load = ref_file_prefix.replace("@",str(i))
            bpfile_prefix = bgen_to_load.replace(".bgen","")
        else:
            bgen_to_load = ref_file_prefix
            bpfile_prefix = bgen_to_load.replace(".bgen","") + f".{i}"

        #check if the file exists
        if convert=="bfile":
            is_file_exist = os.path.exists(bpfile_prefix+".bed")
            make_flag = "--make-bed"
        else:
            is_file_exist = os.path.exists(bpfile_prefix+".pgen")
            make_flag = "--make-pgen vzs"

        #figure out sample file
        if sample is not None:
            if sample=="auto":
                sample_flag = f"--sample {bpfile_prefix}.sample"
            else:
                sample_flag = f"--sample {sample}"
        else:
            sample_flag=""

        # figure memory
        if memory is not None:
            memory_flag = f"--memory {memory}"
        else:
            memory_flag = ""

        #if not existing or overwrite is True
        conversion_success = True
        if (not is_file_exist) or overwrite:
            script_vcf_to_bfile = f"""
            {plink2} \
                --bgen {bgen_to_load} {bgen_mode} {sample_flag}\
                --chr {i} \
                {make_flag} \
                --rm-dup force-first \
                --threads {threads}{memory_flag}\
                --out {bpfile_prefix}
            """
            # execute conversion
            try:
                if convert=="bfile":
                    log.write(f"  -Converting BGEN to bfile: {bpfile_prefix}.bim/bed/fam...")
                else:
                    log.write(f"  -Converting BGEN to pfile: {bpfile_prefix}.pgen/pvar/psam...")
                output = subprocess.check_output(script_vcf_to_bfile, stderr=subprocess.STDOUT, shell=True,text=True)
                plink_log+=output + "\n"
                # Verify conversion succeeded by checking if output file exists
                if convert=="bfile":
                    conversion_success = os.path.exists(bpfile_prefix+".bed")
                else:
                    conversion_success = os.path.exists(bpfile_prefix+".pgen")
                # Read PLINK log file if it exists
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file) as f:
                        plink_log += f.read() + "\n"
            except subprocess.CalledProcessError as e:
                log.write(e.output)
                conversion_success = False
                # Try to read PLINK log file even on error
                plink_log_file = bpfile_prefix + ".log"
                if os.path.exists(plink_log_file):
                    with open(plink_log_file) as f:
                        plink_log += f.read() + "\n"
        else:
            log.write(f"  -PLINK {convert} for CHR {i} exists. Skipping...")

        # Load BIM/PVAR data if requested and conversion was successful
        # Note: load_bim must be True to populate ref_bims, otherwise it remains empty
        if load_bim == True:
            # Only load if file exists (either from successful conversion or pre-existing)
            if convert == "bfile":
                if os.path.exists(bpfile_prefix+".bim"):
                    ref_bims = _load_single_bim_to_ref_bims(bpfile_prefix, ref_bims, log)
                else:
                    log.write(f"  -Warning: BIM file not found for CHR {i}: {bpfile_prefix}.bim. Skipping load.")
            elif os.path.exists(bpfile_prefix+".pvar") or os.path.exists(bpfile_prefix+".pvar.zst"):
                ref_bims = _load_single_pvar_to_ref_bims(bpfile_prefix, ref_bims, log)
            else:
                log.write(f"  -Warning: PVAR file not found for CHR {i}: {bpfile_prefix}.pvar. Skipping load.")
    return ref_file_prefix_converted, plink_log, ref_bims
