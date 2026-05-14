import gzip
from collections.abc import Mapping
from typing import Any, Dict, Optional, Union

import pandas as pd

from gwaslab.bd.bd_common_data import get_format_dict, get_formats_list
from gwaslab.info.g_Log import Log


def _pre_rename_dtype_map(
    meta_data: Mapping[str, Any], dtypes: Mapping[str, Any]
) -> dict[Union[str, int], Any]:
    """Map formatbook ``format_datatype`` keys to labels pandas used before rename."""
    if "format_header" in meta_data:
        fh = meta_data["format_header"]
        if fh is None or fh is False:
            out: dict[Union[str, int], Any] = {}
            for k, v in dtypes.items():
                try:
                    out[int(k)] = v
                except (TypeError, ValueError):
                    out[str(k)] = v
            return out
    return dict(dtypes)


def _count_leading_lines_with_prefix(path: str, prefix: str) -> int:
    """Count initial lines starting with ``prefix`` (e.g. skip VCF/PLINK2 ``##`` meta)."""
    n = 0
    if path.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(prefix):
                    n += 1
                else:
                    break
    else:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(prefix):
                    n += 1
                else:
                    break
    return n


def _multiline_header_skiprows(
    meta_data: Mapping[str, Any],
    load_kwargs_dict: dict[str, Any],
    user_kwargs: Mapping[str, Any],
) -> None:
    """
    Skip trailing header rows after the column-name row (formatbook ``format_header_lines``).

    Line 0 of the file (after any leading ``skiprows``) is the column header; lines
    ``lead + 1`` … ``lead + format_header_lines - 1`` are additional header lines to skip.
    """
    if "skiprows" in user_kwargs:
        return
    raw_nl = meta_data.get("format_header_lines", 1)
    try:
        n_header_lines = int(raw_nl)
    except (TypeError, ValueError):
        return
    if n_header_lines <= 1:
        return
    if "format_header" not in meta_data:
        return
    fh = meta_data["format_header"]
    if fh is None or fh is False:
        return

    def _header_row_index(skiprows_val: Any) -> int | None:
        if skiprows_val is None:
            return 0
        if isinstance(skiprows_val, int):
            return skiprows_val
        if isinstance(skiprows_val, list):
            if not skiprows_val:
                return 0
            if skiprows_val == list(range(len(skiprows_val))):
                return len(skiprows_val)
        return None

    existing = load_kwargs_dict.get("skiprows")
    lead = _header_row_index(existing)
    if lead is None:
        return
    extra = list(range(lead + 1, lead + n_header_lines))
    if not extra:
        return
    if existing is None:
        load_kwargs_dict["skiprows"] = extra
    elif isinstance(existing, int):
        load_kwargs_dict["skiprows"] = list(range(existing)) + extra
    else:
        load_kwargs_dict["skiprows"] = list(existing) + extra


def _read_tabular(path: str, fmt: str, **kwargs: Any) -> pd.DataFrame:

    # default
    load_kwargs_dict = {"sep":"\t",
                      "header":None}

    # if specified by user
    if len(kwargs)>0:
        load_kwargs_dict = kwargs

    # load format
    meta_data, rename_dictionary = get_format_dict(fmt)

    if "format_separator" in meta_data and "sep" not in kwargs:
        load_kwargs_dict["sep"] = meta_data["format_separator"]

    # format_comment: single char -> pandas ``comment=``; multi-char -> *line prefix* at file start:
    # count consecutive leading lines starting with that string and set ``skiprows`` (no ``comment=``),
    # since pandas only allows length-1 ``comment`` and "#" would remove the "#CHROM" header row.
    if "format_comment" in meta_data and meta_data["format_comment"] is not None:
        fc = meta_data["format_comment"]
        if (
            isinstance(fc, str)
            and len(fc) > 1
            and "skiprows" not in kwargs
        ):
            skip_n = _count_leading_lines_with_prefix(path, fc)
            if skip_n:
                load_kwargs_dict["skiprows"] = list(range(skip_n))
        elif isinstance(fc, str) and len(fc) == 1 and "comment" not in kwargs:
            load_kwargs_dict["comment"] = fc

    if "format_header" in meta_data and "header" not in kwargs:
        if meta_data["format_header"] is True:
            load_kwargs_dict["header"] = "infer"
        elif meta_data["format_header"] is False:
            load_kwargs_dict["header"] = None
        else:
            load_kwargs_dict["header"] = meta_data["format_header"]

    if "format_na" in meta_data and "na_values" not in kwargs:
        if  meta_data["format_na"] is not None:
            load_kwargs_dict["na_values"] = meta_data["format_na"]

    _multiline_header_skiprows(meta_data, load_kwargs_dict, kwargs)

    #######################################################################################
    df = pd.read_csv(path, **load_kwargs_dict)
    #######################################################################################

    # format_datatype keys = on-disk names before rename; only cast columns that exist
    if "format_datatype" in meta_data:
        dtype_map = _pre_rename_dtype_map(meta_data, meta_data["format_datatype"])
        dtype_map = {k: v for k, v in dtype_map.items() if k in df.columns}
        if dtype_map:
            df = df.astype(dtype_map)

    # rename columns
    # False or None => no header row was read; pandas columns are 0,1,2,... and format_dict
    # keys are "0","1",... — convert to int so rename matches. True => file supplies names;
    # rename_dictionary keys match those header strings.
    if "format_header" in meta_data:
        fh = meta_data["format_header"]
        if fh is None or fh is False:
            num_to_name = {int(k): v for k, v in rename_dictionary.items()}
            df = df.rename(columns=num_to_name)
        else:
            df = df.rename(columns=rename_dictionary)
    else:
        df = df.rename(columns=rename_dictionary)

    return df

def read_bim(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_bim")
    return df

def read_fam(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_fam")
    return df

def read_psam(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_psam")
    return df

def read_pvar(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="plink_pvar")
    return df

def read_bgen_sample(path: str) -> pd.DataFrame:
    df = _read_tabular(path,fmt="bgen_sample")
    return df
