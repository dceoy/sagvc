#!/usr/bin/env python

import sys
from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta
from luigi.util import requires

from .core import SagvcTask


@requires(FetchReferenceFasta)
class CreateBiallelicSnpVcf(SagvcTask):
    input_vcf_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 90

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        input_vcf = Path(self.input_vcf_path)
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(
                    Path(input_vcf.stem).stem + f'.biallelic_snp.vcf.gz{s}'
                )
            ) for s in ['', '.tbi']
        ]

    def run(self):
        run_id = Path(self.input_vcf_path).stem
        self.print_log(f'Create a common biallelic SNP VCF:\t{run_id}')
        input_vcf = Path(self.input_vcf_path).resource()
        fa = Path(self.fa_path).resource()
        biallelic_snp_vcf = Path(self.output()[0].path)
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=biallelic_snp_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} SelectVariants'
                + f' --variant {input_vcf}'
                + f' --reference {fa}'
                + f' --output {biallelic_snp_vcf}'
                + ' --select-type-to-include SNP'
                + ' --restrict-alleles-to BIALLELIC'
                + ' --lenient'
            ),
            input_files_or_dirs=[input_vcf, fa],
            output_files_or_dirs=[
                biallelic_snp_vcf, f'{biallelic_snp_vcf}.tbi'
            ]
        )


@requires(FetchReferenceFasta)
class CreateWgsIntervalList(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.fa_path).stem + '.wgs.interval_list'
            )
        )

    def run(self):
        fa = Path(self.fa_path).resolve()
        run_id = fa.stem
        self.print_log(f'Create a WGS interval list:\t{run_id}')
        output_interval = Path(self.output().path)
        dest_dir = output_interval.parent
        raw_interval = dest_dir.joinpath(f'{fa.stem}.raw.interval_list')
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=dest_dir, **self.sh_config,
            env={'JAVA_TOOL_OPTIONS': '-Xmx{}m'.format(int(self.memory_mb))}
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} ScatterIntervalsByNs'
                + f' --REFERENCE {fa}'
                + ' --OUTPUT_TYPE ACGT'
                + f' --OUTPUT {raw_interval}'
            ),
            input_files_or_dirs=fa, output_files_or_dirs=raw_interval
        )
        self.run_shell(
            args=(
                'set -e && grep'
                + ' -e \'^@\' -e \'^chr[0-9XYM]\\+\\s\''
                + f' {raw_interval} > {output_interval}'
            ),
            input_files_or_dirs=raw_interval,
            output_files_or_dirs=output_interval
        )
        self.remove_files_and_dirs(raw_interval)


@requires(CreateWgsIntervalList)
class CreateWgsExclusionIntervalListBed(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    samtools = luigi.Parameter(default='samtools')
    gatk = luigi.Parameter(default='gatk')
    bedtools = luigi.Parameter(default='bedtools')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        interval_list = Path(self.input().path)
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{interval_list.stem}.excl.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        yield CreateExclusionIntervalListBed(
            interval_list_path=self.input().path, fa_path=self.fa_path,
            dest_dir_path=self.dest_dir_path, pigz=self.pigz,
            pbzip2=self.pbzip2, samtools=self.samtools, gatk=self.gatk,
            bedtools=self.bedtools, bgzip=self.bgzip, tabix=self.tabix,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            sh_config=self.sh_config
        )


class CreateIntervalListBed(SagvcTask):
    interval_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        interval_list = Path(self.input().path)
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{interval_list.stem}.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        run_id = Path(self.input().path).stem
        self.print_log(f'Create an interval_list BED:\t{run_id}')
        interval_list = Path(self.input().path).resolve()
        bed = Path(self.output()[0].path)
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/interval_list2bed.py'
        )
        self.setup_shell(
            run_id=run_id, commands=[self.bgzip, self.tabix],
            cwd=interval_list.parent, **self.sh_config
        )
        self.run_shell(
            args=(
                'set -eo pipefail'
                + f' && {sys.executable} {pyscript} {interval_list}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {bed}'
            ),
            input_files_or_dirs=interval_list, output_files_or_dirs=bed
        )
        self.tabix_tbi(tsv_path=bed, tabix=self.tabix, preset='bed')


@requires(CreateIntervalListBed, FetchReferenceFasta)
class CreateExclusionIntervalListBed(SagvcTask):
    interval_list_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    bedtools = luigi.Parameter(default='bedtools')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        interval_list = Path(self.input().path)
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{interval_list.stem}.excl.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        input_bed = Path(self.input()[0][0].path)
        run_id = input_bed.stem
        self.print_log(f'Create an exclusion interval_list BED:\t{run_id}')
        fai = Path(f'{self.fa_path}.fai').resolve()
        excl_bed = Path(self.output()[0].path)
        genome_bed_path = self.output()[2].path
        self.setup_shell(
            run_id=run_id, commands=[self.bgzip, self.tabix],
            cwd=input_bed.parent, **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input; '
                    '[print("{0}\\t0\\t{1}".format(*s.split()[:2]))'
                    ' for s in input()];'
                ) + f' {fai}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {genome_bed_path}'
            ),
            input_files_or_dirs=fai, output_files_or_dirs=genome_bed_path
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {self.bedtools} subtract'
                + f' -a {genome_bed_path} -b {input_bed}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {excl_bed}'
            ),
            input_files_or_dirs=[genome_bed_path, input_bed],
            output_files_or_dirs=excl_bed
        )
        self.remove_files_and_dirs(genome_bed_path)
        for p in [genome_bed_path, excl_bed]:
            self.tabix_tbi(tsv_path=p, tabix=self.tabix, preset='bed')


class CreateRegionListBed(SagvcTask):
    region_list_path = luigi.Parameter()
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        region_list = Path(self.region_list_path)
        return [
            luigi.LocalTarget(
                region_list.parent.joinpath(region_list.stem + f'.bed.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        region_list = Path(self.input().path)
        run_id = region_list.stem
        self.print_log(f'Create a region list BED:\t{run_id}')
        bed = Path(self.output()[0].path)
        self.setup_shell(
            run_id=run_id, commands=[self.bgzip, self.tabix],
            cwd=region_list.parent, **self.sh_config
        )
        pycmd = (
            'from fileinput import input;'
            '[(lambda a, b, c: print(f"{a}\t{int(b) - 1}\t{c}"))'
            '(*s.strip().replace(":", "-").split("-")) for s in input()];'
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {sys.executable}'
                + f' -c \'{pycmd}\' {region_list}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {bed}'
            ),
            input_files_or_dirs=region_list, output_files_or_dirs=bed
        )
        self.tabix_tbi(tsv_path=bed, tabix=self.tabix, preset='bed')


class CreateIntervalListWithBed(SagvcTask):
    bed_path = luigi.Parameter()
    seq_dict_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        return luigi.LocalTarget(
            dest_dir.joinpath(Path(self.bed_path).stem + '.interval_list')
        )

    def run(self):
        interval_list = Path(self.output().path)
        run_id = interval_list.stem
        self.print_log(f'Create an interval_list file:\t{run_id}')
        bed = Path(self.bed_path).resolve()
        seq_dict = Path(self.seq_dict_path).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=interval_list.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} BedToIntervalList'
                + f' --INPUT {bed}'
                + f' --SEQUENCE_DICTIONARY {seq_dict}'
                + f' --OUTPUT {interval_list}'
            ),
            input_files_or_dirs=[bed, seq_dict],
            output_files_or_dirs=interval_list
        )


@requires(CreateIntervalListBed)
class UncompressIntervalListBed(SagvcTask):
    bgzip = luigi.Parameter(default='bgzip')
    dest_dir_path = luigi.Parameter(default='.')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.input()[0].path).stem
            )
        )

    def run(self):
        output_bed = Path(self.output().path)
        run_id = output_bed.stem
        self.print_log(f'Uncompress bgzip BED:\t{run_id}')
        bed_gz = Path(self.input()[0].path)
        self.setup_shell(
            run_id=run_id, commands=self.bgzip, cwd=output_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.bgzip} -@ {self.n_cpu} -dc {bed_gz}'
                + f' > {output_bed}'
            ),
            input_files_or_dirs=bed_gz, output_files_or_dirs=output_bed
        )


if __name__ == '__main__':
    luigi.run()
