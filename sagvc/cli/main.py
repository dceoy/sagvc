#!/usr/bin/env python
"""
Somatic and Germline Variant Calling Pipeline

Usage:
    sagvc download [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--use-gnomad-v3]
        [--dest-dir=<path>]
    sagvc write-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--src-path=<path>|--src-url=<url>] [--dest-dir=<path>]
    sagvc haplotypecaller [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--dest-dir=<path>]
        [--interval-list=<interval_list_path>|--bed=<bed_path>] <fa_path>
        <dbsnp_vcf_path> <sam_path>...
    sagvc mutect2 [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--dest-dir=<path>]
        [--interval-list=<interval_list_path>|--bed=<bed_path>]
        [--tumor-sample=<name>] [--normal-sample=<name>] <fa_path>
        <germline_resource_vcf_path> <sam_path>...
    sagvc -h|--help
    sagvc --version

Commands:
    download                Download and preprocess hg38 resources
    haplotypecaller         Identify germline short variants
    mutect2                 Identify somatic short variants

Options:
    -h, --help              Print help and exit
    --version               Print version and exit
    --debug, --info         Execute a command with debug|info messages
    --cpus=<int>            Limit CPU cores used
    --workers=<int>         Specify the maximum number of workers [default: 1]
    --skip-cleaning         Skip incomlete file removal when a task fails
    --print-subprocesses    Print STDOUT/STDERR outputs from subprocesses
    --use-gnomad-v3         Use gnomAD v3 instead of v2
    --dest-dir=<path>       Specify a destination directory path [default: .]
    --src-path=<path>       Specify a source path
    --src-url=<url>         Specify a source URL
    --interval-list=<interval_list_path>, --bed=<bed_path>
                            Specify a path to an interval_list or BED
    --tumor-sample=<name>   Specify a tumor sample name
    --normal-sample=<name>  Specify a normal sample name

Args:
    <fa_path>               Path to an reference FASTA file
                            (The index and sequence dictionary are required.)
    <dbsnp_vcf_path>        Path to a dbSNP VCF file
    <germline_resource_vcf_path>
                            Path to a germline resource VCF file
    <sam_path>              Path to a sorted CRAM or BAM file
"""

import logging
import os
from math import floor
from pathlib import Path

from docopt import docopt
from ftarc.cli.util import (build_luigi_tasks, fetch_executable, print_log,
                            read_yml)
from psutil import cpu_count, virtual_memory

from .. import __version__
from ..task.downloader import (DownloadAndProcessResourceFiles,
                               WritePassingAfOnlyVcf)


def main():
    args = docopt(__doc__, version=__version__)
    if args['--debug']:
        log_level = 'DEBUG'
    elif args['--info']:
        log_level = 'INFO'
    else:
        log_level = 'WARNING'
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=log_level
    )
    logger = logging.getLogger(__name__)
    logger.debug(f'args:{os.linesep}{args}')
    print_log(f'Start the workflow of sagvc {__version__}')
    n_cpu = int(args['--cpus'] or cpu_count())
    memory_mb = virtual_memory().total / 1024 / 1024 / 2
    sh_config = {
        'log_dir_path': args['--dest-dir'],
        'remove_if_failed': (not args['--skip-cleaning']),
        'quiet': (not args['--print-subprocesses']),
        'executable': fetch_executable('bash')
    }
    if args['download']:
        url_dict = read_yml(
            path=Path(__file__).parent.parent.joinpath('static/urls.yml')
        )
        command_dict = {
            'bwa': fetch_executable(
                'bwa-mem2' if args['--use-bwa-mem2'] else 'bwa'
            ),
            'msisensor_pro': fetch_executable('msisensor-pro'),
            **{
                c: fetch_executable(c) for c in [
                    'wget', 'pbzip2', 'bgzip', 'pigz', 'samtools', 'tabix',
                    'gatk', 'bedtools'
                ]
            }
        }
        anns = (
            {k for k in ['snpeff', 'funcotator', 'vep'] if args[f'--{k}']}
            or {'snpeff', 'funcotator', 'vep'}
        )
        n_worker = min(int(args['--workers']), len(anns), n_cpu)
        n_cpu_per_worker = max(1, floor(n_cpu / n_worker))
        memory_mb = int(
            virtual_memory().total / 1024 / 1024 / 2 / n_worker
        )
        common_kwargs = {
            'dest_dir_path': args['--dest-dir'], 'sh_config': sh_config
        }
        build_luigi_tasks(
            tasks=[
                DownloadAndProcessResourceFiles(
                    src_url_dict=url_dict,
                    use_gnomad_exome=args['--use-gnomad-exome'],
                    use_bwa_mem2=args['--use-bwa-mem2'], **command_dict,
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb,
                    **common_kwargs
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['write-af-only-vcf']:
        dest_dir_path = str(Path(args['--dest-dir']).resolve())
        build_luigi_tasks(
            tasks=[
                WritePassingAfOnlyVcf(
                    src_path=(
                        str(Path(args['--src-path']).resolve())
                        if args['--src-path'] else ''
                    ),
                    src_url=(args['--src-url'] or ''),
                    dest_dir_path=dest_dir_path,
                    **{c: fetch_executable(c) for c in ['wget', 'bgzip']},
                    n_cpu=n_cpu, sh_config=sh_config
                )
            ],
            log_level=log_level
        )
    elif args['haplotypecaller']:
        pass
    elif args['mutect2']:
        pass
