#!/usr/bin/env python

import re
import sys
from pathlib import Path

import luigi
from luigi.util import requires

from .core import SagvcTask


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
        output_vcf = Path(self.dest_dir_path).resolve().joinpath(
            Path(Path(self.input_vcf_path).stem).stem + '.biallelic_snp.vcf.gz'
        )
        return [luigi.LocalTarget(f'{output_vcf}{s}') for s in ['', '.tbi']]

    def run(self):
        run_id = Path(self.input_vcf_path).stem
        self.print_log(f'Create a common biallelic SNP VCF:\t{run_id}')
        input_vcf = Path(self.input_vcf_path).resolve()
        fa = Path(self.fa_path).resolve()
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
        output_interval_list = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_interval_list.parent,
            **self.sh_config,
            env={'JAVA_TOOL_OPTIONS': '-Xmx{}m'.format(int(self.memory_mb))}
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} ScatterIntervalsByNs'
                + f' --REFERENCE {fa}'
                + ' --OUTPUT_TYPE ACGT'
                + f' --OUTPUT {output_interval_list}'
            ),
            input_files_or_dirs=fa, output_files_or_dirs=output_interval_list
        )


@requires(CreateWgsIntervalList)
class CreateWgsIntervalListBeds(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    bedtools = luigi.Parameter(default='bedtools')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        interval_list = Path(self.input().path)
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{interval_list.stem}.{s}')
            ) for s in [
                'bed.gz', 'bed.gz.tbi', 'bed', 'excl.bed.gz',
                'excl.bed.gz.tbi', 'excl.bed'
            ]
        ]

    def run(self):
        yield [
            CreateIntervalListBed(
                interval_list_path=self.input().path,
                dest_dir_path=self.dest_dir_path, bgzip=self.bgzip,
                tabix=self.tabix, n_cpu=self.n_cpu, sh_config=self.sh_config
            ),
            CreateExclusionIntervalListBed(
                interval_list_path=self.input().path, fa_path=self.fa_path,
                dest_dir_path=self.dest_dir_path, bedtools=self.bedtools,
                bgzip=self.bgzip, tabix=self.tabix, n_cpu=self.n_cpu,
                sh_config=self.sh_config
            )
        ]


class CreateIntervalListBed(SagvcTask):
    interval_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        bed = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.interval_list_path).stem + '.bed'
        )
        return [luigi.LocalTarget(f'{bed}{s}') for s in ['.gz', '.gz.tbi', '']]

    def run(self):
        bed = Path(self.output()[2].path)
        run_id = bed.stem
        self.print_log(f'Create an interval_list BED:\t{run_id}')
        interval_list = Path(self.interval_list_path).resolve()
        bed_gz = Path(self.output()[0].path)
        pyscript = Path(__file__).resolve().parent.parent.joinpath(
            'script/interval_list2bed.py'
        )
        self.setup_shell(
            run_id=run_id, commands=[self.bgzip, self.tabix],
            cwd=interval_list.parent, **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {sys.executable}'
                + f' {pyscript} {interval_list} > {bed}'
            ),
            input_files_or_dirs=interval_list, output_files_or_dirs=bed
        )
        self.run_shell(
            args=f'set -e && {self.bgzip} -@ {self.n_cpu} -c {bed} > {bed_gz}',
            input_files_or_dirs=bed, output_files_or_dirs=bed_gz
        )
        self.tabix_tbi(tsv_path=bed_gz, tabix=self.tabix, preset='bed')


@requires(CreateIntervalListBed)
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
        excl_bed = Path(self.dest_dir_path).resolve().joinpath(
            Path(Path(self.input()[0].path).stem).stem + '.excl.bed'
        )
        return [
            luigi.LocalTarget(f'{excl_bed}{s}') for s in ['.gz', '.gz.tbi', '']
        ]

    def run(self):
        input_bed = Path(self.input()[0].path)
        run_id = Path(input_bed.stem).stem
        self.print_log(f'Create an exclusion interval_list BED:\t{run_id}')
        fai = Path(f'{self.fa_path}.fai').resolve()
        excl_bed_gz = Path(self.output()[0].path)
        excl_bed = Path(self.output()[2].path)
        dest_dir = excl_bed.parent
        genome_bed = dest_dir.joinpath(f'{fai.stem}.bed')
        self.setup_shell(
            run_id=run_id, commands=[self.bedtools, self.bgzip, self.tabix],
            cwd=dest_dir, **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {sys.executable}'
                + ' -c \'{}\''.format(
                    'from fileinput import input; '
                    '[print("{0}\\t0\\t{1}".format(*s.split()[:2]))'
                    ' for s in input()];'
                ) + f' {fai} > {genome_bed}'
            ),
            input_files_or_dirs=fai, output_files_or_dirs=genome_bed
        )
        self.run_shell(
            args=(
                f'set -eo pipefail && {self.bedtools} subtract'
                + f' -a {genome_bed} -b {input_bed}'
                + f' | {self.bgzip} -@ {self.n_cpu} -c > {excl_bed}'
            ),
            input_files_or_dirs=[genome_bed, input_bed],
            output_files_or_dirs=excl_bed
        )
        self.remove_files_and_dirs(genome_bed)
        self.run_shell(
            args=(
                f'set -e && {self.bgzip} -@ {self.n_cpu} -c {excl_bed}'
                + f' > {excl_bed_gz}'
            ),
            input_files_or_dirs=excl_bed, output_files_or_dirs=excl_bed_gz
        )
        self.tabix_tbi(tsv_path=excl_bed_gz, tabix=self.tabix, preset='bed')


class CreateRegionListBed(SagvcTask):
    region_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        bed = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.region_list_path).stem + '.bed'
        )
        return [luigi.LocalTarget(f'{bed}{s}') for s in ['.gz', '.gz.tbi', '']]

    def run(self):
        run_id = Path(self.region_list_path).stem
        self.print_log(f'Create a region list BED:\t{run_id}')
        region_list = Path(self.region_list_path).resolve()
        bed_gz = Path(self.output()[0].path)
        bed = Path(self.output()[2].path)
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
                + f' -c \'{pycmd}\' {region_list} > {bed}'
            ),
            input_files_or_dirs=region_list, output_files_or_dirs=bed
        )
        self.run_shell(
            args=f'set -e && {self.bgzip} -@ {self.n_cpu} -c {bed} > {bed_gz}',
            input_files_or_dirs=bed, output_files_or_dirs=bed_gz
        )
        self.tabix_tbi(tsv_path=bed_gz, tabix=self.tabix, preset='bed')


class CreateIntervalListWithBed(SagvcTask):
    bed_path = luigi.Parameter()
    fa_dict_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        return luigi.LocalTarget(
            dest_dir.joinpath(
                Path(re.sub(r'\.gz', '', self.bed_path)).stem
                + '.interval_list'
            )
        )

    def run(self):
        interval_list = Path(self.output().path)
        run_id = interval_list.stem
        self.print_log(f'Create an interval_list file:\t{run_id}')
        bed = Path(self.bed_path).resolve()
        fa_dict = Path(self.fa_dict_path).resolve()
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
                + f' --SEQUENCE_DICTIONARY {fa_dict}'
                + f' --OUTPUT {interval_list}'
            ),
            input_files_or_dirs=[bed, fa_dict],
            output_files_or_dirs=interval_list
        )


class SplitEvaluationIntervals(SagvcTask):
    fa_path = luigi.Parameter()
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    scatter_count = luigi.IntParameter(default=1)
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def requires(self):
        if self.interval_list_path:
            return super().requires()
        else:
            return CreateWgsIntervalList(
                fa_path=self.fa_path,
                dest_dir_path=str(
                    Path(self.dest_dir_path).joinpath(
                        Path(self.fa_path).stem + '.wgs'
                    )
                ),
                gatk=self.gatk, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )

    def output(self):
        input_interval_list = Path(
            self.interval_list_path or self.input().path
        ).resolve()
        if self.scatter_count > 1:
            run_dir = Path(self.dest_dir_path).resolve().joinpath(
                f'{input_interval_list.stem}.split_in_{self.scatter_count}'
            )
            return [
                luigi.LocalTarget(
                    run_dir.joinpath(f'{i:04d}-scattered.interval_list')
                ) for i in range(self.scatter_count)
            ]
        else:
            return [luigi.LocalTarget(input_interval_list.resolve())]

    def run(self):
        input_interval_list = Path(
            self.interval_list_path or self.input().path
        ).resolve()
        run_id = input_interval_list.stem
        output_intervals = [Path(o.path) for o in self.output()]
        scatter_count = len(output_intervals)
        self.print_log(f'Split an interval list in {scatter_count}:\t{run_id}')
        fa = Path(self.fa_path).resolve()
        run_dir = output_intervals[0].parent
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} SplitIntervals'
                + f' --reference {fa}'
                + f' --intervals {input_interval_list}'
                + f' --scatter-count {scatter_count}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=[input_interval_list, fa],
            output_files_or_dirs=[*output_intervals, run_dir]
        )


@requires(CreateBiallelicSnpVcf)
class CreateBiallelicSnpIntervalList(SagvcTask):
    fa_path = luigi.Parameter()
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 90

    def output(self):
        snp_vcf = Path(self.input()[0].path)
        return (
            [
                luigi.LocalTarget(
                    snp_vcf.parent.joinpath(f'{snp_vcf.stem}.interval_list')
                )
            ] + self.input()
        )

    def run(self):
        output_interval_list = Path(self.output()[0].path)
        run_id = output_interval_list.stem
        self.print_log(
            f'Create a common biallelic SNP interval_list:\t{run_id}'
        )
        snp_vcf = Path(self.input()[0].path)
        fa = Path(self.fa_path).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_interval_list.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} VcfToIntervalList'
                + f' --INPUT {snp_vcf}'
                + f' --REFERENCE_SEQUENCE {fa}'
                + f' --OUTPUT {output_interval_list}'
            ),
            input_files_or_dirs=[snp_vcf, fa],
            output_files_or_dirs=output_interval_list
        )


class IntersectBed(SagvcTask):
    input_bed_paths = luigi.ListParameter()
    output_bed_path = luigi.Parameter()
    bedtools = luigi.Parameter(default='bedtools')
    sh_config = luigi.DictParameter(default=dict())
    priority = 70

    def output(self):
        return luigi.LocalTarget(
            Path(self.output_bed_path).resolve()
        )

    def run(self):
        assert len(self.input_bed_paths) > 1
        output_bed = Path(self.output().path)
        run_id = Path(output_bed.stem).stem
        self.print_log(f'Create an intersect BED:\t{run_id}')
        input_beds = [Path(p).resolve() for p in self.input_bed_paths]
        self.setup_shell(
            run_id=run_id, commands=self.bedtools, cwd=output_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.bedtools} intersect'
                + f' -a {input_beds[0]} -b'
                + ''.join([f' {b}' for b in input_beds[1:]])
                + f' > {output_bed}'
            ),
            input_files_or_dirs=input_beds, output_files_or_dirs=output_bed
        )


if __name__ == '__main__':
    luigi.run()
