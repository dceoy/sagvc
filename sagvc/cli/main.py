#!/usr/bin/env python
"""
Somatic and Germline Variant Calling Pipeline

Usage:
    sagvc download [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--use-gnomad-v3]
        [--use-msisensor-pro] [--dest-dir=<path>]
    sagvc target [--debug|--info] [--cpus=<int>] [--workers=<int>]
        [--skip-cleaning] [--print-subprocesses] [--dest-dir=<path>]
        [--cnv-blacklist=<path>] <fa_path> <bed_path>
    sagvc write-af-only-vcf [--debug|--info] [--cpus=<int>]
        [--src-path=<path>|--src-url=<url>] [--dest-dir=<path>]
    sagvc haplotypecaller [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>]
        [--interval-list=<path>] [--dbsnp-vcf=<path>] (--resource-vcf=<path>)
        <fa_path> <normal_sam_path>
    sagvc mutect2 [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>]
        [--interval-list=<path>] [--biallelic-snp-vcf=<path>]
        [--germline-resource-vcf=<path>] [--tumor-sample=<name>]
        [--normal-sample=<name>] <fa_path> <tumor_sam_path> <normal_sam_path>
    sagvc delly [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>] [--excl-bed=<path>]
        [--tumor-sample=<name>] [--normal-sample=<name>] <fa_path>
        <tumor_sam_path> <normal_sam_path>
    sagvc manta [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>] [--exome] [--bed=<path>]
        <fa_path> <tumor_sam_path> <normal_sam_path>
    sagvc cnvkit [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>] [--seq-method=<type>]
        [--access-bed=<path>] [--refflat-txt=<path>] <fa_path> <tumor_sam_path>
        <normal_sam_path>
    sagvc callcopyratiosegments [--debug|--info] [--cpus=<int>]
        [--skip-cleaning] [--print-subprocesses] [--dest-dir=<path>]
        [--snp-interval-list=<path>] [--preproc-interval-list=<path>]
        <fa_path> <tumor_sam_path> <normal_sam_path>
    sagvc msisensor [--debug|--info] [--cpus=<int>] [--skip-cleaning]
        [--print-subprocesses] [--dest-dir=<path>] [--use-msisensor-pro]
        [--bed=<path>] [--microsatellites-tsv=<path>] <fa_path>
        <tumor_sam_path> <normal_sam_path>
    sagvc -h|--help
    sagvc --version

Commands:
    download                Download and preprocess hg38 resources
    target                  Prepare targeted region resources
    write-af-only-vcf       Extract and write only AF from VCF INFO
    haplotypecaller         Call germline short variants using GATK
    mutect2                 Call somatic short variants using GATK
    delly                   Call somatic structural variants using Delly
    manta                   Call somatic structural variants using Manta
    cnvkit                  Call somatic CNVs using CNVkit
    callcopyratiosegments   Call somatic CNVs using GATK
    msisensor               Evaluate MSI using MSIsensor

Options:
    -h, --help              Print help and exit
    --version               Print version and exit
    --debug, --info         Execute a command with debug|info messages
    --cpus=<int>            Limit CPU cores used
    --workers=<int>         Specify the maximum number of workers [default: 1]
    --skip-cleaning         Skip incomlete file removal when a task fails
    --print-subprocesses    Print STDOUT/STDERR outputs from subprocesses
    --use-gnomad-v3         Use gnomAD v3 instead of v2
    --use-msisensor-pro     Use MSIsensor-pro instead of MSIsensor
    --dest-dir=<path>       Specify a destination directory path [default: .]
    --cnv-blacklist=<path>  Specify a path to a CNV blacklist file for GATK
    --src-path=<path>       Specify a source file path
    --src-url=<url>         Specify a source URL
    --interval-list=<path>  Specify a path to an interval_list file
    --dbsnp-vcf=<path>      Specify a path to a dbSNP VCF file
    --resource-vcf=<path>   Specify a path to a known SNP and INDEL VCF file
    --biallelic-snp-vcf=<path>
                            Specify a path to a common biallelic SNP VCF file
    --germline-resource-vcf=<path>
                            Specify a path to a germline population VCF file
    --excl-bed=<path>       Specify a path to an exclusion BED file
    --tumor-sample=<name>   Specify a tumor sample name
    --normal-sample=<name>  Specify a normal sample name
    --exome                 Set options for WES input
    --bed=<path>            Specify a path to a BED file
    --seq-method=<type>     Specify a sequencing assay type [default: wgs]
    --access-bed=<path>     Specify a path to a CNV accessible region BED file
    --refflat-txt=<path>    Specify a path to a refFlat text file
    --snp-interval-list=<path>
                            Specify a path to a common sites VCF file
    --preproc-interval-list=<path>
                            Specify a path to a preprocessed interval_list file
    --microsatellites-tsv=<path>
                            Specify a path to a microsatellites TSV file

Args:
    <fa_path>               Path to a reference FASTA file
                            (The index and sequence dictionary are required.)
    <bed_path>              Path to a sorted target BED file
    <tumor_sam_path>        Path to a sorted tumor CRAM file
    <normal_sam_path>       Path to a sorted normal CRAM file
"""

import logging
import os
from math import ceil, floor
from pathlib import Path

from docopt import docopt
from ftarc.cli.util import (build_luigi_tasks, fetch_executable, print_log,
                            print_yml, read_yml)
from psutil import cpu_count, virtual_memory

from .. import __version__
from ..task.callcopyratiosegments import CallCopyRatioSegmentsTumor
from ..task.cnvkit import CallSomaticCnvWithCnvkit
from ..task.delly import CallSomaticStructualVariantsWithDelly
from ..task.downloader import (DownloadAndProcessResourceFiles,
                               PrepareTargetedReosourceFiles,
                               WritePassingAfOnlyVcf)
from ..task.haplotypecaller import FilterVariantTranches
from ..task.manta import CallSomaticStructualVariantsWithManta
from ..task.msisensor import ScoreMsiWithMsisensor
from ..task.mutect2 import FilterMutectCalls


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
    if args['download'] or args['target']:
        n_worker = int(args['--workers'])
    elif args['haplotypecaller'] or args['mutect2']:
        n_worker = n_cpu
    elif args['callcopyratiosegments']:
        n_worker = min(n_cpu, 2)
    else:
        n_worker = 1
    n_cpu_per_worker = max(floor(n_cpu / n_worker), 1)
    memory_mb_per_worker = ceil(
        virtual_memory().total / 1024 / 1024 / 2 / n_worker
    )
    print_yml([
        {'n_worker': n_worker}, {'n_cpu_per_worker': n_cpu_per_worker},
        {'memory_mb_per_worker': memory_mb_per_worker}
    ])
    sh_config = {
        'log_dir_path': args['--dest-dir'],
        'remove_if_failed': (not args['--skip-cleaning']),
        'quiet': (not args['--print-subprocesses']),
        'executable': fetch_executable('bash')
    }
    if args['download']:
        build_luigi_tasks(
            tasks=[
                DownloadAndProcessResourceFiles(
                    src_url_dict=read_yml(
                        path=Path(__file__).parent.parent.joinpath(
                            'static/urls.yml'
                        )
                    ),
                    dest_dir_path=args['--dest-dir'],
                    **{
                        c: fetch_executable(c) for c in [
                            'wget', 'pbzip2', 'bgzip', 'pigz', 'samtools',
                            'tabix', 'gatk', 'bedtools'
                        ]
                    },
                    cnvkitpy=fetch_executable('cnvkit.py'),
                    msisensor=fetch_executable(
                        'msisensor-pro' if args['--use-msisensor-pro']
                        else 'msisensor'
                    ),
                    use_gnomad_v3=args['--use-gnomad-v3'],
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['target']:
        build_luigi_tasks(
            tasks=[
                PrepareTargetedReosourceFiles(
                    bed_path=args['<bed_path>'], fa_path=args['<fa_path>'],
                    cnv_blacklist_path=(args['--cnv-blacklist'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    **{
                        c: fetch_executable(c)
                        for c in ['gatk', 'bgzip', 'tabix',  'bedtools']
                    },
                    cnvkitpy=fetch_executable('cnvkit.py'),
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['write-af-only-vcf']:
        build_luigi_tasks(
            tasks=[
                WritePassingAfOnlyVcf(
                    src_path=(
                        str(Path(args['--src-path']).resolve())
                        if args['--src-path'] else ''
                    ),
                    src_url=(args['--src-url'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    **{c: fetch_executable(c) for c in ['wget', 'bgzip']},
                    n_cpu=n_cpu_per_worker, sh_config=sh_config
                )
            ],
            log_level=log_level
        )
    elif args['haplotypecaller']:
        assert bool(args['--resource-vcf']), '--resource-vcf required'
        build_luigi_tasks(
            tasks=[
                FilterVariantTranches(
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    dbsnp_vcf_path=(args['--dbsnp-vcf'] or ''),
                    resource_vcf_paths=args['--resource-vcf'],
                    interval_list_path=(args['--interval-list'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    gatk=fetch_executable('gatk'),
                    samtools=fetch_executable('samtools'),
                    save_memory=(memory_mb_per_worker < 8192),
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config, scatter_count=n_worker
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['mutect2']:
        for k in ['--biallelic-snp-vcf', '--germline-resource-vcf']:
            assert bool(args[k]), f'{k} required'
        build_luigi_tasks(
            tasks=[
                FilterMutectCalls(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    tumor_sample_name=args['--tumor-sample'],
                    normal_sample_name=args['--normal-sample'],
                    biallelic_snp_vcf_path=args['--biallelic-snp-vcf'],
                    germline_resource_vcf_path=args['--germline-resource-vcf'],
                    interval_list_path=(args['--interval-list'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    gatk=fetch_executable('gatk'),
                    samtools=fetch_executable('samtools'),
                    save_memory=(memory_mb_per_worker < 8192),
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config, scatter_count=n_worker
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['delly']:
        build_luigi_tasks(
            tasks=[
                CallSomaticStructualVariantsWithDelly(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    tumor_sample_name=args['--tumor-sample'],
                    normal_sample_name=args['--normal-sample'],
                    excl_bed_path=(args['--excl-bed'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    delly=fetch_executable('delly'),
                    bcftools=fetch_executable('bcftools'),
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['manta']:
        build_luigi_tasks(
            tasks=[
                CallSomaticStructualVariantsWithManta(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'], bed_path=(args['--bed'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    python2=fetch_executable('python2'),
                    configmantapy_path=fetch_executable('configManta.py'),
                    exome=args['--exome'], n_cpu=n_cpu_per_worker,
                    memory_mb=memory_mb_per_worker, sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['cnvkit']:
        build_luigi_tasks(
            tasks=[
                CallSomaticCnvWithCnvkit(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    access_bed_path=(args['--access-bed'] or ''),
                    refflat_txt_path=(args['--refflat-txt'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    cnvkitpy=fetch_executable('cnvkit.py'),
                    samtools=fetch_executable('samtools'),
                    rscript=fetch_executable('Rscript'),
                    seq_method=args['--seq-method'],
                    n_cpu=n_cpu_per_worker, sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['callcopyratiosegments']:
        for k in ['--snp-interval-list', '--preproc-interval-list']:
            assert bool(args[k]), f'{k} required'
        build_luigi_tasks(
            tasks=[
                CallCopyRatioSegmentsTumor(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    snp_interval_list_path=args['--snp-interval-list'],
                    preproc_interval_list_path=args['--preproc-interval-list'],
                    dest_dir_path=args['--dest-dir'],
                    gatk=fetch_executable('gatk'),
                    save_memory=(memory_mb_per_worker < 8192),
                    n_cpu=n_cpu_per_worker, memory_mb=memory_mb_per_worker,
                    sh_config=sh_config,
                )
            ],
            workers=n_worker, log_level=log_level
        )
    elif args['msisensor']:
        build_luigi_tasks(
            tasks=[
                ScoreMsiWithMsisensor(
                    tumor_cram_path=args['<tumor_sam_path>'],
                    normal_cram_path=args['<normal_sam_path>'],
                    fa_path=args['<fa_path>'],
                    microsatellites_tsv_path=(
                        args['--microsatellites-tsv'] or ''
                    ),
                    bed_path=(args['--bed'] or ''),
                    dest_dir_path=args['--dest-dir'],
                    msisensor=fetch_executable(
                        'msisensor-pro' if args['--use-msisensor-pro']
                        else 'msisensor'
                    ),
                    samtools=fetch_executable('samtools'),
                    n_cpu=n_cpu_per_worker, sh_config=sh_config
                )
            ],
            workers=n_worker, log_level=log_level
        )
