"""
GWAS Catalog summary statistics downloader.

Workflow overview
-----------------
This module downloads summary statistics from GWAS Catalog by GCST accession
ID (for example, ``GCST90270926``). The workflow is intentionally defensive
because different studies expose files via slightly different routes.

1) Validate input and prepare output
   - Validate GCST format (must match ``GCST`` + digits).
   - If ``output_dir`` is not provided, use GWASLab default data directory
     (``options.paths["data_directory"]``), then create directory if needed.
   - Create a study-specific subdirectory ``<output_dir>/<GCST_ID>/`` so all
     related files are stored together.
   - Delay final filename selection until a concrete download URL is chosen.

2) Construct FTP bucket path from GCST number
   - Convert the numeric part of GCST ID to the 1,000-study bucket used by EBI.
   - Example: ``GCST90270926`` maps to
     ``GCST90270001-GCST90271000/GCST90270926/``.

3) Discover downloadable files from FTP index pages (preferred path)
   - Read the study directory listing.
   - If ``harmonised=True`` and ``harmonised/`` exists, inspect it first.
   - Select the first likely sumstats file (``.tsv.gz``, ``.txt.gz``, ``.vcf.gz``,
     and uncompressed variants), excluding metadata/checksum files.
   - If no harmonised file is found, fall back to non-harmonised files in the
     study root directory. If ``harmonised=False``, select from raw files only.

4) Fallback discovery routes
   - If FTP listing is unavailable or does not provide a usable file, try direct
     summary-statistics API download endpoints.
   - If still unresolved, query GWAS metadata endpoints and recursively extract
     URL-like fields to find a sumstats file link.

5) Download sidecar files and verify integrity when possible
   - If available, download ``<main_file>-meta.yaml`` and ``md5sum.txt``.
   - If ``md5sum.txt`` includes checksum for the selected file, verify local
     MD5 after download.

6) Normalize URL and download safely
   - Normalize ``ftp://`` and ``http://ftp.ebi.ac.uk`` links to HTTPS.
   - Infer output filename from the selected URL when ``filename`` is not set.
   - Stream file in chunks to ``*.part`` and atomically rename on success.
   - Respect ``overwrite=False`` by returning existing output path unchanged.

7) Register downloaded data in local records
   - Add/update a ``downloaded`` config record keyed by GCST ID.
   - Include main path and sidecar/checksum metadata when available.

This priority order ensures harmonised files are preferred whenever available,
while still providing robust fallback behavior across GWAS Catalog endpoints.
"""

import hashlib
import os
import re
from typing import Any, List, Optional
from urllib.parse import urljoin, urlparse

import requests

from gwaslab.bd.bd_config import options
from gwaslab.bd.bd_download import update_description
from gwaslab.info.g_Log import Log


def _looks_like_sumstats_url(url: str, gcst_id: str) -> bool:
    """Check whether a URL likely points to a summary statistics resource."""
    if not isinstance(url, str):
        return False
    lowered = url.lower()
    gcst_lower = gcst_id.lower()
    return (
        (lowered.startswith("http://") or lowered.startswith("https://") or lowered.startswith("ftp://"))
        and gcst_lower in lowered
        and (
            "summary_statistics" in lowered
            or "summary-statistics" in lowered
            or "sumstats" in lowered
            or lowered.endswith(".tsv.gz")
            or lowered.endswith(".txt.gz")
            or lowered.endswith(".vcf.gz")
        )
    )


def _collect_urls(payload: Any) -> list[str]:
    """Recursively collect URL-like strings from nested JSON objects."""
    urls: list[str] = []
    if isinstance(payload, dict):
        for value in payload.values():
            urls.extend(_collect_urls(value))
    elif isinstance(payload, list):
        for value in payload:
            urls.extend(_collect_urls(value))
    elif isinstance(payload, str):
        if payload.startswith("http://") or payload.startswith("https://") or payload.startswith("ftp://"):
            urls.append(payload)
    return urls


def _normalize_download_url(url: str) -> str:
    """Convert FTP URLs to HTTPS URLs so requests can download them."""
    if url.startswith("ftp://ftp.ebi.ac.uk/"):
        return "https://ftp.ebi.ac.uk/" + url.split("ftp://ftp.ebi.ac.uk/", 1)[1]
    if url.startswith("http://ftp.ebi.ac.uk/"):
        return "https://ftp.ebi.ac.uk/" + url.split("http://ftp.ebi.ac.uk/", 1)[1]
    return url


def _build_ftp_gcst_dir(gcst_id: str) -> str:
    """Build GWAS Catalog FTP directory URL for a GCST accession."""
    gcst_numeric = gcst_id.replace("GCST", "")
    gcst_number = int(gcst_numeric)
    width = max(len(gcst_numeric), 8)
    block_start = ((gcst_number - 1) // 1000) * 1000 + 1
    block_end = block_start + 999
    block_label = "GCST{start:0{w}d}-GCST{end:0{w}d}".format(
        start=block_start,
        end=block_end,
        w=width,
    )
    return (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        f"{block_label}/{gcst_id}/"
    )


def _list_directory_entries(url: str, timeout: int) -> list[str]:
    """List href targets from an EBI FTP directory index page."""
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    hrefs = re.findall(r'href="([^"]+)"', response.text)
    return [h for h in hrefs if h not in ("../",)]


def _list_directory_entries_safe(url: str, timeout: int) -> list[str]:
    """List href targets and return an empty list on request errors."""
    try:
        return _list_directory_entries(url, timeout=timeout)
    except requests.exceptions.RequestException:
        return []


def _is_likely_sumstats_file(name: str) -> bool:
    """Check whether filename looks like a downloadable sumstats file."""
    lowered = name.lower()
    if lowered.endswith("/"):
        return False
    if lowered.endswith(("-meta.yaml", ".md5", "md5sum.txt", ".sha256", ".tbi", ".log")):
        return False
    return lowered.endswith((".tsv.gz", ".txt.gz", ".vcf.gz", ".tsv", ".txt", ".vcf"))


def _pick_best_sumstats_entry(entries: list[str], gcst_id: str, prefer_harmonised: bool) -> str | None:
    """Pick the best sumstats filename from directory entries."""
    candidates = [entry for entry in entries if _is_likely_sumstats_file(entry)]
    if not candidates:
        return None

    gcst_lower = gcst_id.lower()
    lowered = [entry.lower() for entry in candidates]

    if prefer_harmonised:
        preferred_patterns = [
            f"{gcst_lower}.h.tsv.gz",
            f"{gcst_lower}.h.txt.gz",
            f"{gcst_lower}.h.vcf.gz",
            f"{gcst_lower}.h.tsv",
            f"{gcst_lower}.h.txt",
            f"{gcst_lower}.h.vcf",
        ]
        for pattern in preferred_patterns:
            for idx, name in enumerate(lowered):
                if name == pattern:
                    return candidates[idx]
        for idx, name in enumerate(lowered):
            if ".h." in name:
                return candidates[idx]
    else:
        preferred_patterns = [
            f"{gcst_lower}.tsv.gz",
            f"{gcst_lower}.txt.gz",
            f"{gcst_lower}.vcf.gz",
            f"{gcst_lower}.tsv",
            f"{gcst_lower}.txt",
            f"{gcst_lower}.vcf",
        ]
        for pattern in preferred_patterns:
            for idx, name in enumerate(lowered):
                if name == pattern:
                    return candidates[idx]

    return candidates[0]


def _compute_md5(path: str, chunk_size: int = 1024 * 1024) -> str:
    """Compute MD5 checksum for a local file."""
    md5 = hashlib.md5()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            md5.update(chunk)
    return md5.hexdigest()


def _parse_md5sum_for_file(md5_text: str, target_name: str) -> str | None:
    """Parse md5sum.txt content and return checksum for target filename."""
    target_name = target_name.strip()
    lines = [line.strip() for line in md5_text.splitlines() if line.strip()]
    if not lines:
        return None

    parsed_pairs = []
    for line in lines:
        # Typical md5sum format: "<md5>  <filename>" or "<md5> *<filename>"
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        checksum = parts[0].strip().lower()
        filename = parts[1].strip().lstrip("*")
        filename = filename.replace("\\", "/").split("/")[-1]
        parsed_pairs.append((filename, checksum))

    for filename, checksum in parsed_pairs:
        if filename == target_name:
            return checksum

    # Fallback: if checksum file has exactly one parsed row, use it.
    if len(parsed_pairs) == 1:
        return parsed_pairs[0][1]
    return None


def _download_if_exists(
    base_url: str,
    entry_name: str,
    output_dir: str,
    timeout: int,
    chunk_size: int,
    overwrite: bool,
    log: Log,
    verbose: bool
) -> str | None:
    """Download a sidecar file from directory index if present."""
    entries = _list_directory_entries_safe(base_url, timeout=timeout)
    if entry_name not in entries:
        return None

    target_url = _normalize_download_url(urljoin(base_url, entry_name))
    target_path = os.path.abspath(os.path.join(output_dir, entry_name))
    if os.path.exists(target_path) and not overwrite:
        log.write(f" -Auxiliary file exists: {target_path} (overwrite=False)", verbose=verbose)
        return target_path

    log.write(f" -Downloading auxiliary file: {target_url}", verbose=verbose)
    _stream_download(target_url, target_path, timeout=timeout, chunk_size=chunk_size)
    return target_path


def _prepare_gcst_output_dir(base_output_dir: str, gcst_id: str, log: Log, verbose: bool) -> str:
    """Prepare a GCST-specific output directory, handling name collisions safely."""
    gcst_output_dir = os.path.abspath(os.path.join(base_output_dir, gcst_id))
    if os.path.exists(gcst_output_dir) and not os.path.isdir(gcst_output_dir):
        fallback_dir = gcst_output_dir + "_files"
        log.warning(
            f" -Path exists as a file, using fallback directory instead: {fallback_dir}",
            verbose=verbose,
        )
        gcst_output_dir = fallback_dir

    os.makedirs(gcst_output_dir, exist_ok=True)
    return gcst_output_dir


def _discover_ftp_sumstats_url(
    gcst_id: str,
    timeout: int,
    log: Log,
    verbose: bool,
    harmonised: bool = True
) -> str | None:
    """Discover sumstats URL from FTP index, with optional harmonised preference."""
    study_dir_url = _build_ftp_gcst_dir(gcst_id)
    try:
        log.write(f" -Trying FTP directory: {study_dir_url}", verbose=verbose)
        root_entries = _list_directory_entries(study_dir_url, timeout=timeout)
    except requests.exceptions.RequestException:
        return None

    # 1) Optionally prioritize harmonised directory.
    if harmonised:
        harmonised_path = None
        for entry in root_entries:
            if entry.rstrip("/").lower() == "harmonised":
                harmonised_path = entry
                break

        if harmonised_path is not None:
            harmonised_url = urljoin(study_dir_url, harmonised_path)
            try:
                harm_entries = _list_directory_entries(harmonised_url, timeout=timeout)
                picked = _pick_best_sumstats_entry(harm_entries, gcst_id=gcst_id, prefer_harmonised=True)
                if picked is not None:
                    return _normalize_download_url(urljoin(harmonised_url, picked))
            except requests.exceptions.RequestException:
                pass

    # 2) Raw files in the study directory (or fallback when harmonised missing).
    picked = _pick_best_sumstats_entry(root_entries, gcst_id=gcst_id, prefer_harmonised=False)
    if picked is not None:
        return _normalize_download_url(urljoin(study_dir_url, picked))

    return None


def _discover_sumstats_url(gcst_id: str, timeout: int, log: Log, verbose: bool) -> str | None:
    """Try several metadata endpoints and extract a downloadable sumstats URL."""
    metadata_endpoints = [
        f"https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/{gcst_id}",
        f"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/{gcst_id}",
        f"https://www.ebi.ac.uk/gwas/rest/api/studies/{gcst_id}",
    ]

    headers = {"Accept": "application/json"}

    for endpoint in metadata_endpoints:
        try:
            log.write(f" -Trying metadata endpoint: {endpoint}", verbose=verbose)
            response = requests.get(endpoint, headers=headers, timeout=timeout)
            if response.status_code != 200:
                continue

            data = response.json()
            for url in _collect_urls(data):
                if _looks_like_sumstats_url(url, gcst_id):
                    return _normalize_download_url(url)
        except requests.exceptions.RequestException:
            continue
        except ValueError:
            continue

    return None


def _stream_download(url: str, output_path: str, timeout: int, chunk_size: int) -> None:
    """Download URL content to output_path with streaming."""
    tmp_path = output_path + ".part"
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with open(tmp_path, "wb") as handle:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    handle.write(chunk)
    os.replace(tmp_path, output_path)


def download_sumstats(
    gcst_id: str,
    output_dir: str | None = None,
    filename: str | None = None,
    harmonised: bool = True,
    timeout: int = 60,
    chunk_size: int = 1024 * 1024,
    overwrite: bool = False,
    verbose: bool = True,
    log: Log | None = None,
) -> str:
    """
    Download GWAS Catalog summary statistics for a given GCST accession.

    Parameters
    ----------
    gcst_id : str
        GWAS Catalog study accession (e.g. "GCST90002446").
    output_dir : str, optional
        Directory to store downloaded file. If None, uses the GWASLab default
        data directory from ``options.paths["data_directory"]``. Files are then
        stored under ``<output_dir>/<GCST_ID>/``.
    filename : str, optional
        Output filename. If None, defaults to "{GCST}_sumstats.tsv.gz".
    harmonised : bool, optional
        If True, prioritize harmonised files. If harmonised files are not
        available, fall back to raw files. If False, select raw files directly
        (default: True).
    timeout : int, optional
        HTTP request timeout in seconds (default: 60).
    chunk_size : int, optional
        Download chunk size in bytes (default: 1MB).
    overwrite : bool, optional
        Whether to overwrite existing output file (default: False).
    verbose : bool, optional
        Whether to print log messages (default: True).
    log : Log, optional
        Logger object. If None, a new logger is created.

    Returns
    -------
    str
        Absolute path of downloaded file.

    Raises
    ------
    ValueError
        If GCST ID format is invalid.
    RuntimeError
        If no downloadable summary statistics URL can be found.
    requests.exceptions.RequestException
        If the download request fails.
    """
    logger = log if log is not None else Log()

    gcst_id = gcst_id.strip().upper()
    if not re.match(r"^GCST\d+$", gcst_id):
        raise ValueError(f"Invalid GCST ID: {gcst_id}. Expected format like GCST90002446.")

    if output_dir is None:
        output_dir = options.paths["data_directory"]
    gcst_output_dir = _prepare_gcst_output_dir(output_dir, gcst_id, log=logger, verbose=verbose)

    logger.write(f" -Searching GWAS Catalog summary statistics for {gcst_id}", verbose=verbose)
    selected_url: str | None = None

    # Prefer FTP-organized study directory and optionally harmonised files.
    selected_url = _discover_ftp_sumstats_url(
        gcst_id,
        timeout=timeout,
        log=logger,
        verbose=verbose,
        harmonised=harmonised,
    )

    direct_download_endpoints = [
        f"https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/{gcst_id}/download",
        f"https://www.ebi.ac.uk/gwas/summary-statistics/api/studies/{gcst_id}/downloads",
    ]

    last_error: Exception | None = None
    if selected_url is None:
        for endpoint in direct_download_endpoints:
            try:
                logger.write(f" -Trying direct download endpoint: {endpoint}", verbose=verbose)
                # Probe endpoint first so we can still pick filename from URL if needed.
                response = requests.get(endpoint, stream=True, timeout=timeout)
                response.raise_for_status()
                selected_url = response.url
                response.close()
                break
            except requests.exceptions.RequestException as exc:
                last_error = exc
                continue

    if selected_url is None:
        selected_url = _discover_sumstats_url(gcst_id, timeout=timeout, log=logger, verbose=verbose)

    if selected_url is None:
        if last_error is not None:
            raise RuntimeError(
                f"Unable to locate downloadable summary statistics for {gcst_id}. "
                f"Tried direct endpoints and metadata discovery; last error: {last_error!s}"
            )
        raise RuntimeError(
            f"Unable to locate downloadable summary statistics for {gcst_id}. "
            "Tried direct endpoints and metadata discovery."
        )

    selected_url = _normalize_download_url(selected_url)
    inferred_name = os.path.basename(urlparse(selected_url).path.rstrip("/"))
    if not inferred_name:
        inferred_name = f"{gcst_id}_sumstats.tsv.gz"

    resolved_filename = filename if filename is not None else inferred_name
    output_path = os.path.abspath(os.path.join(gcst_output_dir, resolved_filename))
    if os.path.exists(output_path) and not overwrite:
        logger.write(f" -File already exists: {output_path} (overwrite=False)", verbose=verbose)
    else:
        logger.write(f" -Using discovered download URL: {selected_url}", verbose=verbose)
        _stream_download(selected_url, output_path, timeout=timeout, chunk_size=chunk_size)
        logger.write(f" -Downloaded summary statistics to {output_path}", verbose=verbose)

    # Download optional sidecar files and perform MD5 verification when possible.
    selected_dir_url = selected_url.rsplit("/", 1)[0] + "/"
    candidate_dirs = [selected_dir_url]
    selected_dir_path = urlparse(selected_dir_url).path.rstrip("/")
    if selected_dir_path.endswith("/harmonised"):
        parent_url = selected_dir_url.rstrip("/").rsplit("/", 1)[0] + "/"
        candidate_dirs.append(parent_url)

    yaml_downloaded = False
    yaml_path: str | None = None
    md5_path: str | None = None
    for base_dir in candidate_dirs:
        yaml_name = f"{inferred_name}-meta.yaml"
        yaml_path = _download_if_exists(
            base_url=base_dir,
            entry_name=yaml_name,
            output_dir=gcst_output_dir,
            timeout=timeout,
            chunk_size=chunk_size,
            overwrite=overwrite,
            log=logger,
            verbose=verbose,
        )
        if yaml_path is not None:
            yaml_downloaded = True

        md5_candidate = _download_if_exists(
            base_url=base_dir,
            entry_name="md5sum.txt",
            output_dir=gcst_output_dir,
            timeout=timeout,
            chunk_size=chunk_size,
            overwrite=overwrite,
            log=logger,
            verbose=verbose,
        )
        if md5_candidate is not None:
            md5_path = md5_candidate
            break

    if not yaml_downloaded:
        logger.write(" -No metadata YAML found next to selected file", verbose=verbose)

    expected_md5: str | None = None
    md5_verified: bool | None = None
    if md5_path is not None:
        with open(md5_path, encoding="utf-8") as handle:
            md5_text = handle.read()
        expected_md5 = _parse_md5sum_for_file(md5_text, inferred_name)
        if expected_md5 is not None:
            observed_md5 = _compute_md5(output_path, chunk_size=chunk_size)
            if observed_md5.lower() != expected_md5.lower():
                md5_verified = False
                raise RuntimeError(
                    f"MD5 mismatch for {os.path.basename(output_path)}: expected {expected_md5}, observed {observed_md5}"
                )
            md5_verified = True
            logger.write(f" -MD5 check passed for {os.path.basename(output_path)}", verbose=verbose)
        else:
            logger.write(f" -md5sum.txt found but no checksum entry for {inferred_name}", verbose=verbose)
    else:
        logger.write(" -No md5sum.txt found for checksum validation", verbose=verbose)

    # Register downloaded dataset in config, following bd_download conventions.
    try:
        file_format = "unknown"
        lower_path = output_path.lower()
        if lower_path.endswith(".tsv.gz"):
            file_format = "tsv.gz"
        elif lower_path.endswith(".txt.gz"):
            file_format = "txt.gz"
        elif lower_path.endswith(".vcf.gz"):
            file_format = "vcf.gz"
        else:
            ext = os.path.splitext(output_path)[1].lstrip(".")
            if ext:
                file_format = ext

        record = {
            "local_path": output_path,
            "description": f"GWAS Catalog summary statistics for {gcst_id}",
            "suggested_use": "Use as external GWAS summary statistics for downstream analyses.",
            "source": "GWAS Catalog",
            "url": selected_url,
            "format": file_format,
            "gcst_id": gcst_id,
            "harmonised_requested": harmonised,
            "harmonised_downloaded": ".h." in inferred_name.lower(),
        }
        if yaml_path is not None:
            record["meta_yaml"] = {"local_path": yaml_path}
        if md5_path is not None:
            record["md5_file"] = {"local_path": md5_path}
        if expected_md5 is not None:
            record["md5sum"] = expected_md5
        if md5_verified is not None:
            record["md5_verified"] = md5_verified

        update_description(gcst_id, record, log=logger)
    except Exception as e:
        logger.warning(
            f" -Downloaded file but failed to update local record: {e!s}",
            verbose=verbose,
        )

    return output_path

