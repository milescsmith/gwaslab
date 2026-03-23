"""
GWASLab CLI - Main entry point.

Simple flat CLI structure with flags for operations.
"""

import argparse
import sys
from typing import Optional

from gwaslab_cli.handlers import (
    load_sumstats,
    run_config,
    run_config_show,
    run_path,
    run_fb_list,
    run_fb_show,
    run_fb_update,
    run_version,
    run_download,
)


def main(argv: Optional[list] = None) -> None:
    """Main entry point for GWASLab CLI."""
    if argv is None:
        argv = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="gwaslab",
        description=(
            "GWASLab: A Python package for processing GWAS summary statistics\n"
            "\n"
            "Processing order (when multiple flags are given):\n"
            "  1. QC            --qc / --remove / --remove-dup / --normalize\n"
            "  2. Harmonize     --harmonize [--ref-seq ...]\n"
            "  3. Assign rsID   --assign-rsid  (auto fix_chr+fix_pos if QC not run)\n"
            "  4. rsID→CHR:POS  --rsid-to-chrpos\n"
            "  5. Infer build   --infer-build\n"
            "  6. Liftover      --liftover FROM TO  (auto fix_chr+fix_pos if QC not run)\n"
            "  7. Plot          --plot TYPE         (auto fix_chr+fix_pos if QC not run; writes to --output, then exits)\n"
            "  8. Extract       --extract TYPE      (writes to --output, then exits)\n"
            "  9. Save          --output FILE [--to-fmt FORMAT]"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process sumstats
  gwaslab --input sumstats.tsv --qc --output cleaned.tsv
  gwaslab --input sumstats.tsv --harmonize --ref-seq ref.fa --output harmonized.tsv
  gwaslab --input sumstats.tsv --liftover 19 38 --output lifted_hg38.tsv
  gwaslab --input sumstats.tsv --infer-build --liftover 19 38 --output lifted_hg38.tsv
  gwaslab --input sumstats.tsv --output converted.ldsc --to-fmt ldsc

  # Plot
  gwaslab --input sumstats.tsv --plot manhattan --output manhattan.png
  gwaslab --input sumstats.tsv --plot qq --output qq.png

  # Extract
  gwaslab --input sumstats.tsv --extract lead --output leads.tsv
  gwaslab --input sumstats.tsv --extract novel --efo EFO_0004330 --output novel.tsv

  # Other commands
  gwaslab version
  gwaslab config
  gwaslab config show config
  gwaslab config show reference
  gwaslab path config
  gwaslab formatbook list
  gwaslab download-sumstats GCST90270926
        """,
    )

    # Input/Output
    parser.add_argument("--input", "-i", help="Input sumstats file path")
    parser.add_argument("--out", "--output", "-o", dest="output", help="Output file path")
    parser.add_argument("--fmt", "-f", default="auto", help="Input format (default: auto)")
    parser.add_argument("--to-fmt", default="gwaslab", help="Output format (default: gwaslab)")
    parser.add_argument("--nrows", type=int, help="Number of rows to read (for testing)")

    # Processing flags
    parser.add_argument("--qc", action="store_true", help="Perform quality control")
    parser.add_argument("--harmonize", action="store_true", help="Perform harmonization")
    parser.add_argument(
        "--assign-rsid",
        action="store_true",
        help="Assign rsID to variants (runs fix_chr+fix_pos if basic_check was not run)",
    )
    parser.add_argument("--rsid-to-chrpos", action="store_true", help="Convert rsID to CHR:POS")
    parser.add_argument(
        "--infer-build",
        action="store_true",
        help="Infer genome build from HapMap3 SNP coordinates (hg19/hg38)",
    )
    parser.add_argument(
        "--liftover",
        nargs=2,
        metavar=("FROM_BUILD", "TO_BUILD"),
        help="Liftover from one build to another (runs fix_chr+fix_pos if basic_check was not run)",
    )

    # Plot flag
    parser.add_argument(
        "--plot",
        choices=["manhattan", "qq", "mqq", "regional", "miami", "forest"],
        help="Generate a plot (runs fix_chr+fix_pos if basic_check was not run)",
    )

    # Extract flag
    parser.add_argument("--extract", choices=["lead", "novel", "proxy"], help="Extract variants")

    # Common options
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress output")

    # QC options
    parser.add_argument("--remove", action="store_true", help="Remove bad quality variants")
    parser.add_argument("--remove-dup", action="store_true", help="Remove duplicated variants")
    parser.add_argument("--normalize", action="store_true", help="Normalize indels")

    # Harmonization options
    parser.add_argument("--ref-seq", help="Reference FASTA file")
    parser.add_argument("--ref-rsid-tsv", help="Reference rsID HDF5 file")
    parser.add_argument("--ref-rsid-vcf", help="Reference rsID VCF file")
    parser.add_argument("--ref-infer", help="Reference VCF for strand inference")
    parser.add_argument("--ref-maf-threshold", type=float, default=0.4)
    parser.add_argument("--maf-threshold", type=float, default=0.40)
    parser.add_argument("--sweep-mode", action="store_true")
    parser.add_argument("--basic-check", action="store_true", help="Enable basic check during harmonization")
    parser.add_argument(
        "--no-basic-check", action="store_true", help="Disable basic check during harmonization"
    )

    # Output options
    parser.add_argument("--no-gzip", action="store_true")
    parser.add_argument("--bgzip", action="store_true")
    parser.add_argument("--tabix", action="store_true")
    parser.add_argument("--hapmap3", action="store_true")
    parser.add_argument("--exclude-hla", action="store_true")
    parser.add_argument("--hla-lower", type=int, help="Lower bound for HLA exclusion (in Mb)")
    parser.add_argument("--hla-upper", type=int, help="Upper bound for HLA exclusion (in Mb)")
    parser.add_argument("--build", default="19")
    parser.add_argument("--chr-prefix", default="")
    parser.add_argument("--xymt-number", action="store_true", help="Use numeric encoding for XYMT chromosomes")
    parser.add_argument("--n", type=int, help="Sample size to add to output")
    parser.add_argument("--tab-fmt", default="tsv", help="Output table format (tsv, csv, parquet)")
    parser.add_argument("--overwrite", default="empty", help="Overwrite mode for assign-rsid (empty, all, etc.)")

    # Plot options
    parser.add_argument("--sig-level", type=float, default=5e-8)
    parser.add_argument("--ylim", type=float, nargs=2)
    parser.add_argument("--highlight", nargs="+")
    parser.add_argument("--chr", type=int, help="Chromosome for regional plot")
    parser.add_argument("--start", type=int, help="Start position for regional plot")
    parser.add_argument("--end", type=int, help="End position for regional plot")

    # Extract options
    parser.add_argument("--sig-level-extract", type=float, default=5e-8, dest="sig_level_extract")
    parser.add_argument("--windowsizekb", type=int, default=500)
    parser.add_argument("--efo", nargs="+", help="EFO IDs for novel extraction")
    parser.add_argument("--only-novel", action="store_true")

    # Subcommands (version, config, etc.)
    subparsers = parser.add_subparsers(dest="command")

    # Version
    version_parser = subparsers.add_parser("version", help="Show version")
    version_parser.set_defaults(func=run_version)

    # Config
    config_parser = subparsers.add_parser("config", help="Show configuration")
    config_parser.add_argument("--json", action="store_true")
    config_parser.set_defaults(func=run_config)
    config_sub = config_parser.add_subparsers(dest="config_action")
    config_show = config_sub.add_parser("show", help="Show one configured path by keyword")
    config_show.add_argument("keyword", help="Path keyword (e.g., config, reference)")
    config_show.set_defaults(func=run_config_show)

    # Path
    path_parser = subparsers.add_parser("path", help="Resolve path")
    path_parser.add_argument("keyword")
    path_parser.set_defaults(func=run_path)

    # Formatbook
    fb_parser = subparsers.add_parser("formatbook", help="Manage formats")
    fb_sub = fb_parser.add_subparsers(dest="fb_action")
    fb_list = fb_sub.add_parser("list", help="List formats")
    fb_list.add_argument("--json", action="store_true")
    fb_list.set_defaults(func=run_fb_list)
    fb_show = fb_sub.add_parser("show", help="Show format")
    fb_show.add_argument("format")
    fb_show.set_defaults(func=run_fb_show)
    fb_update = fb_sub.add_parser("update", help="Update formatbook")
    fb_update.set_defaults(func=run_fb_update)

    # Download
    dl_parser = subparsers.add_parser("download-sumstats", help="Download from GWAS Catalog")
    dl_parser.add_argument("gcst_id")
    dl_parser.add_argument("--output-dir", "-o")
    dl_parser.set_defaults(func=run_download)

    args = parser.parse_args(argv)

    # Handle subcommands
    if args.command:
        if hasattr(args, "func"):
            args.func(args)
        return

    # Require --input for main operations
    if not args.input:
        parser.error("--input is required (or use a subcommand like 'version', 'config')")

    # Derive input build from --build or --liftover FROM build.
    input_build = args.build
    if args.liftover and args.build == "19" and not args.infer_build:
        input_build = args.liftover[0]

    # Load sumstats
    s = load_sumstats(args.input, args.fmt, args.nrows, build=input_build)
    ran_basic_check = False
    coords_ready = False

    def ensure_coords_ready() -> None:
        """Ensure CHR/POS are standardized for operations that require them."""
        nonlocal coords_ready
        if not coords_ready:
            s.fix_chr(verbose=not args.quiet)
            s.fix_pos(verbose=not args.quiet)
            coords_ready = True

    # Processing
    if args.qc or args.remove or args.remove_dup or args.normalize:
        normalize = args.normalize if args.normalize else (True if args.qc else False)
        s.basic_check(
            remove=args.remove,
            remove_dup=args.remove_dup,
            threads=args.threads,
            normalize=normalize,
            verbose=not args.quiet,
        )
        ran_basic_check = True
        coords_ready = True

    if args.harmonize or args.ref_seq:
        if args.basic_check:
            basic_check = True
        elif args.no_basic_check:
            basic_check = False
        else:
            basic_check = not args.qc
        s.harmonize(
            basic_check=basic_check,
            ref_seq=args.ref_seq,
            ref_rsid_tsv=args.ref_rsid_tsv,
            ref_rsid_vcf=args.ref_rsid_vcf,
            ref_infer=args.ref_infer,
            ref_maf_threshold=args.ref_maf_threshold,
            maf_threshold=args.maf_threshold,
            threads=args.threads,
            sweep_mode=args.sweep_mode,
            verbose=not args.quiet,
        )
        ran_basic_check = ran_basic_check or basic_check
        coords_ready = coords_ready or basic_check

    if args.assign_rsid:
        ensure_coords_ready()
        s.assign_rsid(
            ref_rsid_tsv=args.ref_rsid_tsv,
            ref_rsid_vcf=args.ref_rsid_vcf,
            overwrite=args.overwrite,
            threads=args.threads,
            verbose=not args.quiet,
        )

    if args.rsid_to_chrpos:
        s.rsid_to_chrpos(
            ref_rsid_to_chrpos_vcf=args.ref_rsid_vcf,
            ref_rsid_to_chrpos_hdf5=args.ref_rsid_tsv,
            build=args.build,
            threads=args.threads if args.threads > 1 else 4,
            verbose=not args.quiet,
        )

    if args.infer_build:
        s.infer_build(verbose=not args.quiet)

    if args.liftover:
        ensure_coords_ready()
        from_build, to_build = args.liftover
        s.liftover(
            from_build=from_build,
            to_build=to_build,
            verbose=not args.quiet,
        )

    # Plotting
    if args.plot:
        if not ran_basic_check:
            # Plotting needs sanitized CHR/POS. If basic_check was not run,
            # apply the lightweight fallback of fix_chr + fix_pos.
            ensure_coords_ready()
        if args.plot == "manhattan":
            s.plot_mqq(mode="m", save=args.output, verbose=not args.quiet)
        elif args.plot == "qq":
            s.plot_qq(save=args.output, verbose=not args.quiet)
        elif args.plot == "mqq":
            s.plot_mqq(mode="mqq", save=args.output, verbose=not args.quiet)
        elif args.plot == "regional":
            if not all([args.chr, args.start, args.end]):
                parser.error("--chr, --start, --end required for regional plot")
            s.plot_mqq(mode="r", region=(args.chr, args.start, args.end), save=args.output, verbose=not args.quiet)
        elif args.plot == "miami":
            parser.error("Miami plot requires two inputs. Use Python API: gl.plot_miami2()")
        elif args.plot == "forest":
            s.plot_forest(save=args.output, verbose=not args.quiet)
        return

    # Extraction
    if args.extract:
        if args.extract == "lead":
            result = s.get_lead(
                sig_level=args.sig_level_extract,
                windowsizekb=args.windowsizekb,
                verbose=not args.quiet,
            )
            if args.output:
                result.to_csv(args.output, sep="\t", index=False)
        elif args.extract == "novel":
            result = s.get_novel(
                efo=args.efo or False,
                only_novel=args.only_novel,
                sig_level=args.sig_level_extract,
                verbose=not args.quiet,
            )
            if args.output:
                result.to_csv(args.output, sep="\t", index=False)
        elif args.extract == "proxy":
            parser.error("Proxy extraction not yet implemented in CLI")
        return

    # Output (if not plot/extract which handle their own output)
    if args.output:
        to_format_kwargs = {
            "fmt": args.to_fmt,
            "tab_fmt": args.tab_fmt,
            "gzip": not args.no_gzip,
            "bgzip": args.bgzip,
            "tabix": args.tabix,
            "hapmap3": args.hapmap3,
            "exclude_hla": args.exclude_hla,
            "chr_prefix": args.chr_prefix,
            "xymt_number": args.xymt_number,
            "n": args.n,
            "verbose": not args.quiet,
        }

        if args.hla_lower is not None or args.hla_upper is not None:
            to_format_kwargs["hla_range"] = (args.hla_lower or 25, args.hla_upper or 34)

        s.to_format(args.output, **to_format_kwargs)


if __name__ == "__main__":
    main()

