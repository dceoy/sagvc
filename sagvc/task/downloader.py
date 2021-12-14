#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi
from ftarc.task.downloader import (DownloadAndIndexResourceVcfs,
                                   DownloadResourceFiles)
from ftarc.task.picard import CreateSequenceDictionary
from ftarc.task.samtools import SamtoolsFaidx
from luigi.util import requires

from .callcopyratiosegments import PreprocessIntervals
from .cnvkit import CreateCnvAcccessBed
from .core import SagvcTask
from .msisensor import ScanMicrosatellites
from .resource import (CreateBiallelicSnpIntervalList,
                       CreateExclusionIntervalListBed, CreateIntervalListBed,
                       CreateIntervalListWithBed, CreateWgsIntervalListBeds,
                       IntersectBed)


class PrepareTargetedReosourceFiles(luigi.Task):
    bed_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    cnv_blacklist_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    bedtools = luigi.Parameter(default='bedtools')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 80

    def requires(self):
        fa = Path(self.fa_path).resolve()
        return CreateIntervalListWithBed(
            bed_path=self.bed_path,
            fa_dict_path=str(fa.parent.joinpath(f'{fa.stem}.dict')),
            dest_dir_path=self.dest_dir_path, gatk=self.gatk, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, sh_config=self.sh_config
        )

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        bed_stem = Path(re.sub(r'\.gz', '', self.bed_path)).stem
        return [
            luigi.LocalTarget(dest_dir.joinpath(n)) for n in [
                f'{bed_stem}.interval_list', f'{bed_stem}.bed.gz',
                f'{bed_stem}.bed.gz.tbi', f'{bed_stem}.bed',
                f'{bed_stem}.excl.bed.gz', f'{bed_stem}.excl.bed.gz.tbi',
                f'{bed_stem}.excl.bed', f'{bed_stem}.access.bed',
                *[
                    (
                        bed_stem + (
                            ('.excl.' + Path(self.cnv_blacklist_path).stem)
                            if self.cnv_blacklist_path else ''
                        ) + f'.w{r}s.preproc.interval_list'
                    ) for r in ['x', 'g']
                ]
            ]
        ]

    def run(self):
        interval_list_path = self.input().path
        input_targets = yield [
            CreateIntervalListBed(
                interval_list_path=interval_list_path,
                dest_dir_path=self.dest_dir_path, bgzip=self.bgzip,
                tabix=self.tabix, n_cpu=self.n_cpu, sh_config=self.sh_config
            ),
            CreateExclusionIntervalListBed(
                interval_list_path=interval_list_path, fa_path=self.fa_path,
                dest_dir_path=self.dest_dir_path, bedtools=self.bedtools,
                bgzip=self.bgzip, tabix=self.tabix, n_cpu=self.n_cpu,
                sh_config=self.sh_config
            ),
            CreateCnvAcccessBed(
                fa_path=self.fa_path, dest_dir_path=self.dest_dir_path,
                cnvkitpy=self.cnvkitpy, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            )
        ]
        yield [
            IntersectBed(
                input_bed_paths=[
                    input_targets[2].path, input_targets[0][0].path
                ],
                output_bed_path=self.output()[7].path, bedtools=self.bedtools,
                sh_config=self.sh_config
            ),
            *[
                PreprocessIntervals(
                    fa_path=self.fa_path,
                    cnv_blacklist_path=self.cnv_blacklist_path,
                    interval_list_path=interval_list_path,
                    dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                    exome=bool(i), n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                    sh_config=self.sh_config
                ) for i in range(2)
            ]
        ]


class DownloadAndProcessRegionFiles(luigi.Task):
    src_url_dict = luigi.DictParameter()
    dest_dir_path = luigi.Parameter(default='.')
    wget = luigi.Parameter(default='wget')
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    samtools = luigi.Parameter(default='samtools')
    gatk = luigi.Parameter(default='gatk')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    bedtools = luigi.Parameter(default='bedtools')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    msisensor = luigi.Parameter(default='msisensor')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 80

    def requires(self):
        return [
            DownloadAndIndexReferenceFasta(
                src_url_dict=self.src_url_dict,
                dest_dir_path=self.dest_dir_path, wget=self.wget,
                pigz=self.pigz, pbzip2=self.pbzip2, bgzip=self.bgzip,
                samtools=self.samtools, gatk=self.gatk, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ),
            DownloadAndIndexResourceVcfs(
                src_urls=[
                    self.src_url_dict['cnv_blacklist'],
                    *[
                        v for k, v in self.src_url_dict.items()
                        if k not in {'ref_fa', 'ref_fa_alt', 'cnv_blacklist'}
                    ]
                ],
                dest_dir_path=self.dest_dir_path, run_id='others',
                wget=self.wget, pigz=self.pigz, pbzip2=self.pbzip2,
                bgzip=self.bgzip, tabix=self.tabix, n_cpu=self.n_cpu,
                sh_config=self.sh_config
            )
        ]

    def output(self):
        fa = Path(self.input()[0][0].path)
        cnv_blacklist = Path(self.input()[1][0].path)
        return (
            self.input()[0] + self.input()[1] + [
                luigi.LocalTarget(fa.parent.joinpath(n)) for n in [
                    f'{fa.name}.fai', f'{fa.stem}.dict',
                    f'{fa.stem}.wgs.interval_list', f'{fa.stem}.wgs.bed.gz',
                    f'{fa.stem}.wgs.bed.gz.tbi', f'{fa.stem}.wgs.bed',
                    f'{fa.stem}.wgs.excl.bed.gz',
                    f'{fa.stem}.wgs.excl.bed.gz.tbi',
                    f'{fa.stem}.wgs.excl.bed', f'{fa.stem}.access.bed',
                    *[
                        (
                            f'{fa.stem}.excl.{cnv_blacklist.stem}.w{r}s'
                            + '.preproc.interval_list'
                        ) for r in ['x', 'g']
                    ],
                    f'{fa.stem}.microsatellites.tsv'
                ]
            ]
        )

    def run(self):
        fa_path = self.input()[0][0].path
        dest_dir_path = str(Path(fa_path).parent)
        cnv_blacklist_path = self.input()[1][0].path
        yield [
            CreateWgsIntervalListBeds(
                fa_path=fa_path, dest_dir_path=dest_dir_path, gatk=self.gatk,
                bedtools=self.bedtools, bgzip=self.bgzip, tabix=self.tabix,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            ),
            CreateCnvAcccessBed(
                fa_path=fa_path, dest_dir_path=dest_dir_path,
                cnvkitpy=self.cnvkitpy, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ),
            *[
                PreprocessIntervals(
                    fa_path=fa_path, cnv_blacklist_path=cnv_blacklist_path,
                    dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                    exome=bool(i), n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                    sh_config=self.sh_config
                ) for i in range(2)
            ],
            ScanMicrosatellites(
                fa_path=fa_path, dest_dir_path=dest_dir_path,
                msisensor=self.msisensor, sh_config=self.sh_config
            )
        ]


class DownloadAndIndexReferenceFasta(luigi.Task):
    src_url_dict = luigi.DictParameter()
    dest_dir_path = luigi.Parameter(default='.')
    wget = luigi.Parameter(default='wget')
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    bgzip = luigi.Parameter(default='bgzip')
    samtools = luigi.Parameter(default='samtools')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def requires(self):
        return DownloadResourceFiles(
            src_urls=[
                self.src_url_dict['ref_fa'],
                self.src_url_dict['ref_fa_alt']
            ],
            dest_dir_path=self.dest_dir_path,
            run_id=Path(self.src_url_dict['ref_fa']).stem, wget=self.wget,
            pigz=self.pigz, pbzip2=self.pbzip2, bgzip=self.bgzip,
            n_cpu=self.n_cpu, sh_config=self.sh_config
        )

    def output(self):
        fa = Path(self.input()[0].path)
        return [
            *self.input(), luigi.LocalTarget(f'{fa}.fai'),
            luigi.LocalTarget(fa.parent.joinpath(f'{fa.stem}.dict'))
        ]

    def run(self):
        fa_path = self.input()[0].path
        yield [
            SamtoolsFaidx(
                fa_path=fa_path, samtools=self.samtools,
                sh_config=self.sh_config
            ),
            CreateSequenceDictionary(
                fa_path=fa_path, gatk=self.gatk, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            )
        ]


class DownloadGnomadVcfsAndExtractAf(SagvcTask):
    dest_dir_path = luigi.Parameter(default='.')
    use_gnomad_v3 = luigi.BoolParameter(default=False)
    cloud_storage = luigi.Parameter(default='amazon')
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 100

    def output(self):
        output_vcf = Path(self.dest_dir_path).resolve().joinpath(
            'gnomad.genomes.v3.1.sites.af-only.vcf.gz'
            if self.use_gnomad_v3 else
            'gnomad.exomes.r2.1.1.sites.liftover_grch38.af-only.vcf.gz'
        )
        return [luigi.LocalTarget(f'{output_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_id = Path(Path(output_vcf.stem).stem).stem
        self.print_log(f'Download and process a gnomAD VCF file:\t{run_id}')
        dest_dir = output_vcf.parent
        url_root = {
            'google': 'storage.googleapis.com/gcp-public-data--gnomad',
            'amazon': 'gnomad-public-us-east-1.s3.amazonaws.com',
            'microsoft': 'azureopendatastorage.blob.core.windows.net/gnomad'
        }[self.cloud_storage.lower()]
        if self.use_gnomad_v3:
            urls = [
                (
                    f'https://{url_root}/release/3.1/vcf/genomes/'
                    + f'gnomad.genomes.v3.1.sites.chr{i}.vcf.bgz'
                ) for i in [*range(1, 23), 'X', 'Y']
            ]
        else:
            urls = [
                f'https://{url_root}/release/2.1.1/liftover_grch38/vcf/exomes/'
                + 'gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz'
            ]
        vcf_dict = {
            u: dest_dir.joinpath(Path(Path(u).stem).stem + '.af-only.vcf.gz')
            for u in urls
        }
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/extract_af_only_vcf.py'
        )
        self.setup_shell(
            run_id=run_id,
            commands=[self.wget, self.bgzip, sys.executable, self.gatk],
            cwd=dest_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        for u, v in vcf_dict.items():
            self.run_shell(
                args=(
                    f'set -e && {self.wget} -qSL {u} -O -'
                    + f' | {self.bgzip} -@ {self.n_cpu} -dc'
                    + f' | {sys.executable} {pyscript} -'
                    + f' | {self.bgzip} -@ {self.n_cpu} -c > {v}'
                ),
                output_files_or_dirs=v
            )
        if output_vcf.is_file():
            self.tabix_tbi(tsv_path=output_vcf, tabix=self.tabix, preset='vcf')
        else:
            self.picard_mergevcfs(
                input_vcf_paths=vcf_dict.values(), output_vcf_path=output_vcf,
                picard=self.gatk, remove_input=True
            )


@requires(DownloadGnomadVcfsAndExtractAf, DownloadAndIndexReferenceFasta)
class DownloadAndProcessGnomadVcf(luigi.Task):
    dest_dir_path = luigi.Parameter(default='.')
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    samtools = luigi.Parameter(default='samtools')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 90

    def output(self):
        ba_snp_vcf = Path(self.dest_dir_path).resolve().joinpath(
            Path(Path(self.input()[0][0].path).stem).stem
            + '.biallelic_snp.vcf.gz'
        )
        return (
            self.input()[0] + self.input()[1]
            + [
                luigi.LocalTarget(ba_snp_vcf),
                luigi.LocalTarget(f'{ba_snp_vcf}.tbi'),
                luigi.LocalTarget(
                    ba_snp_vcf.parent.joinpath(
                        f'{ba_snp_vcf.stem}.interval_list'
                    )
                )
            ]
        )

    def run(self):
        yield CreateBiallelicSnpIntervalList(
            input_vcf_path=self.input()[0][0].path,
            fa_path=self.input()[1][0].path, dest_dir_path=self.dest_dir_path,
            gatk=self.gatk, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


@requires(DownloadAndProcessGnomadVcf,  DownloadAndProcessRegionFiles)
class DownloadAndProcessResourceFiles(luigi.WrapperTask):
    priority = 10

    def output(self):
        return self.input()


class WritePassingAfOnlyVcf(SagvcTask):
    src_path = luigi.Parameter(default='')
    src_url = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    wget = luigi.Parameter(default='wget')
    bgzip = luigi.Parameter(default='bgzip')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(Path(self.src_path or self.src_url).stem).stem
                + '.af-only.vcf.gz'
            )
        )

    def run(self):
        assert bool(self.src_path or self.src_url)
        output_vcf = Path(self.output().path)
        run_id = Path(Path(output_vcf.stem).stem).stem
        message = (
            'Write a passing AF-only VCF' if self.src_path
            else 'Download a VCF file and extract passing AF-only records'
        )
        self.print_log(f'{message}:\t{run_id}')
        dest_dir = output_vcf.parent
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/extract_af_only_vcf.py'
        )
        self.setup_shell(
            run_id=run_id,
            commands=[
                *(list() if self.src_path else [self.wget]), self.bgzip,
                sys.executable
            ],
            cwd=dest_dir, **self.sh_config
        )
        if self.src_path:
            src_vcf = Path(self.src_path).resolve()
        else:
            src_vcf = dest_dir.joinpath(Path(self.src_url).name)
            self.run_shell(
                args=f'set -e && {self.wget} -qSL {self.src_url} -O {src_vcf}',
                output_files_or_dirs=src_vcf
            )
        self.run_shell(
            args=(
                f'set -e && {self.bgzip}'
                + f' -@ {self.n_cpu} -dc {src_vcf}'
                + f' | {sys.executable} {pyscript} -'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {output_vcf}'
            ),
            input_files_or_dirs=src_vcf, output_files_or_dirs=output_vcf
        )


if __name__ == '__main__':
    luigi.run()
