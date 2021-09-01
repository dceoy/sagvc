#!/usr/bin/env python

import re
import sys
from pathlib import Path
from socket import gethostname

import luigi
from ftarc.task.downloader import DownloadResourceFiles
from ftarc.task.resource import FetchResourceVcf

from .core import SagvcTask
from .resource import (CreateBiallelicSnpVcf, CreateRegionListBed,
                       CreateWgsIntervalListBeds)


class DownloadAndProcessResourceFiles(luigi.Task):
    src_urls = luigi.ListParameter()
    dest_dir_path = luigi.Parameter(default='.')
    run_id = luigi.Parameter(default=gethostname())
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    bwa = luigi.Parameter(default='bwa')
    samtools = luigi.Parameter(default='samtools')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    gatk = luigi.Parameter(default='gatk')
    bedtools = luigi.Parameter(default='bedtools')
    msisensorpro = luigi.Parameter(default='msisensor-pro')
    use_gnomad_v3 = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    use_bwa_mem2 = luigi.BoolParameter(default=False)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def requires(self):
        return [
            DownloadGnomadVcfsAndExtractAf(
                dest_dir_path=self.dest_dir_path,
                use_gnomad_v3=self.use_gnomad_v3, wget=self.wget,
                bgzip=self.bgzip, gatk=self.gatk, tabix=self.tabix,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            ),
            DownloadResourceFiles(
                src_urls=self.src_urls, dest_dir_path=self.dest_dir_path,
                wget=self.wget, bgzip=self.bgzip, pbzip2=self.pbzip2,
                pigz=self.pigz, n_cpu=self.n_cpu, sh_config=self.sh_config
            )
        ]

    def output(self):
        for i in self.input()[0]:
            yield luigi.LocalTarget(i.path)
            yield luigi.LocalTarget(
                re.sub(r'(vcf\.gz|vcf\.gz\.tbi)$', r'biallelic_snp.\1', i.path)
            )
        for i in self.input()[1]:
            f = Path(i.path)
            yield luigi.LocalTarget(f)
            if f.name.endswith(('.fa', '.fna', '.fasta')):
                for n in [f'{f.name}.fai', f'{f.stem}.dict',
                          f'{f.stem}.wgs.interval_list',
                          f'{f.stem}.wgs.bed.gz',
                          f'{f.stem}.wgs.bed.gz.tbi']:
                    yield luigi.LocalTarget(f.parent.joinpath(n))
            elif f.name.endswith('.vcf.gz'):
                yield luigi.LocalTarget(f'{f}.tbi')

    def run(self):
        gnomad_vcf = Path(self.input()['gnomad_vcf'][0].path)
        downloaded_file_paths = [i.path for i in self.input()[1]]
        fa = [
            Path(p) for p in downloaded_file_paths
            if p.endswith(('.fa', '.fna', '.fasta'))
        ][0]
        yield [
            CreateBiallelicSnpVcf(
                input_vcf_path=str(gnomad_vcf), fa_path=str(fa),
                dest_dir_path=self.dest_dir_path, pigz=self.pigz,
                pbzip2=self.pbzip2, samtools=self.samtools, gatk=self.gatk,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            ),
            CreateWgsIntervalListBeds(
                fa_path=str(fa), dest_dir_path=self.dest_dir_path,
                pigz=self.pigz, pbzip2=self.pbzip2, samtools=self.samtools,
                gatk=self.gatk, bedtools=self.bedtools, bgzip=self.bgzip,
                tabix=self.tabix, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            ),
            *[
                FetchResourceVcf(
                    src_path=p, bgzip=self.bgzip, tabix=self.tabix,
                    n_cpu=self.n_cpu, sh_config=self.sh_config
                ) for p in downloaded_file_paths if p.endswith('.vcf.gz')
            ],
            *[
                CreateRegionListBed(
                    region_list_path=p, bgzip=self.bgzip, tabix=self.tabix,
                    n_cpu=self.n_cpu, sh_config=self.sh_config
                ) for p in downloaded_file_paths if p.endswith('.list')
            ]
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
    priority = 10

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
