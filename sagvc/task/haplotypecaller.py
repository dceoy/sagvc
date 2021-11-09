#!/usr/bin/env python

import re
from itertools import chain
from pathlib import Path

import luigi
from luigi.util import requires

from .core import SagvcTask
from .resource import SplitEvaluationIntervals


@requires(SplitEvaluationIntervals)
class CallVariantsWithHaplotypeCaller(SagvcTask):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dbsnp_vcf_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    samtools = luigi.Parameter(default='samtools')
    add_haplotypecaller_args = luigi.ListParameter(
        default=[
            '--standard-min-confidence-threshold-for-calling', '0',
            '--annotation', 'Coverage', '--annotation', 'ChromosomeCounts',
            '--annotation', 'BaseQuality', '--annotation', 'FragmentLength',
            '--annotation', 'MappingQuality', '--annotation', 'ReadPosition',
            '--create-output-bam-index', 'false'
        ]
    )
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.normal_cram_path).stem
        )
        cram_stem = Path(self.normal_cram_path).stem
        return [
            luigi.LocalTarget(
                run_dir.joinpath(f'{cram_stem}.haplotypecaller.{s}')
            ) for s in ['vcf.gz', 'vcf.gz.tbi', 'cram', 'cram.crai']
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        interval_lists = [Path(i.path) for i in self.input()]
        skip_interval_list_split = (len(interval_lists) == 1)
        fa = Path(self.fa_path).resolve()
        output_path_prefix = '.'.join(str(output_vcf).split('.')[:-2])
        if skip_interval_list_split:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, o.stem)
                for o in interval_lists
            ]
        input_targets = yield [
            HaplotypeCaller(
                input_cram_path=self.normal_cram_path, fa_path=str(fa),
                dbsnp_vcf_path=self.dbsnp_vcf_path,
                interval_list_path=str(o), output_path_prefix=s,
                gatk=self.gatk,
                add_haplotypecaller_args=self.add_haplotypecaller_args,
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ) for o, s in zip(interval_lists, tmp_prefixes)
        ]
        run_id = '.'.join(output_vcf.name.split('.')[:-3])
        self.print_log(
            f'Call germline variants with HaplotypeCaller:\t{run_id}'
        )
        output_cram = Path(self.output()[2].path)
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        if skip_interval_list_split:
            tmp_bam = Path(f'{tmp_prefixes[0]}.bam')
            self.samtools_view(
                input_sam_path=tmp_bam, fa_path=fa,
                output_sam_path=output_cram, samtools=self.samtools,
                n_cpu=self.n_cpu, index_sam=True, remove_input=True
            )
        else:
            tmp_vcfs = [Path(f'{s}.vcf.gz') for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} MergeVcfs'
                    + ''.join(f' --INPUT {v}' for v in tmp_vcfs)
                    + f' --REFERENCE_SEQUENCE {fa}'
                    + f' --OUTPUT {output_vcf}'
                ),
                input_files_or_dirs=[*tmp_vcfs, fa],
                output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
            )
            self.samtools_merge(
                input_sam_paths=[f'{s}.bam' for s in tmp_prefixes],
                fa_path=fa, output_sam_path=output_cram,
                samtools=self.samtools, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, index_sam=True, remove_input=False
            )
            self.remove_files_and_dirs(
                *chain.from_iterable(
                    [o.path for o in t] for t in input_targets
                )
            )


class HaplotypeCaller(SagvcTask):
    input_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    interval_list_path = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    dbsnp_vcf_path = luigi.Parameter(default='')
    gatk = luigi.Parameter(default='gatk')
    add_haplotypecaller_args = luigi.ListParameter(default=list())
    save_memory = luigi.BoolParameter(default=False)
    message = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.{s}')
            for s in ['vcf.gz', 'vcf.gz.tbi', 'bam']
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        input_cram = Path(self.input_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        dbsnp_vcf = (
            Path(self.dbsnp_vcf_path).resolve()
            if self.dbsnp_vcf_path else None
        )
        interval_list = Path(self.interval_list_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        output_vcf = output_files[0]
        run_dir = output_vcf.parent
        self.setup_shell(
            run_id='.'.join(output_vcf.name.split('.')[:-2]),
            commands=self.gatk, cwd=run_dir, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} HaplotypeCaller'
                + f' --input {input_cram}'
                + f' --read-index {input_cram}.crai'
                + f' --reference {fa}'
                + (f' --dbsnp {dbsnp_vcf}' if dbsnp_vcf else '')
                + f' --intervals {interval_list}'
                + f' --native-pair-hmm-threads {self.n_cpu}'
                + ''.join(
                    f' {a}' for a in [
                        *self.add_haplotypecaller_args,
                        *(
                            ['--disable-bam-index-caching', 'true']
                            if self.save_memory else list()
                        )
                    ]
                )
                + f' --output {output_vcf}'
                + f' --bam-output {output_files[2]}'
            ),
            input_files_or_dirs=[
                input_cram, f'{input_cram}.crai', fa, interval_list,
                *([dbsnp_vcf] if dbsnp_vcf else list())
            ],
            output_files_or_dirs=[*output_files, run_dir]
        )


@requires(CallVariantsWithHaplotypeCaller, SplitEvaluationIntervals)
class ScoreVariantsWithCnn(SagvcTask):
    fa_path = luigi.Parameter()
    gatk = luigi.Parameter(default='gatk')
    python = luigi.Parameter(default='python')
    add_cnnscorevariants_args = luigi.ListParameter(
        default=['--tensor-type', 'read_tensor']
    )
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_path_prefix = re.sub(r'\.vcf\.gz$', '', self.input()[0][0].path)
        return [
            luigi.LocalTarget(f'{output_path_prefix}.cnn.vcf.gz{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        interval_lists = [Path(i.path) for i in self.input()[1]]
        skip_interval_list_split = (len(interval_lists) == 1)
        output_vcf = Path(self.output()[0].path)
        output_path_prefix = '.'.join(str(output_vcf).split('.')[:-2])
        if skip_interval_list_split:
            tmp_prefixes = [output_path_prefix]
        else:
            tmp_prefixes = [
                '{0}.{1}'.format(output_path_prefix, o.stem)
                for o in interval_lists
            ]
        input_targets = yield [
            CNNScoreVariants(
                input_vcf_path=self.input()[0][0].path,
                input_cram_path=self.input()[0][2].path, fa_path=self.fa_path,
                interval_list_path=str(o), output_path_prefix=s,
                gatk=self.gatk, python=self.python,
                add_cnnscorevariants_args=self.add_cnnscorevariants_args,
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ) for o, s in zip(interval_lists, tmp_prefixes)
        ]
        run_id = '.'.join(output_vcf.name.split('.')[:-2])
        self.print_log(f'Score variants with CNN:\t{run_id}')
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        if not skip_interval_list_split:
            tmp_vcfs = [Path(f'{s}.vcf.gz') for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} MergeVcfs'
                    + ''.join(f' --INPUT {v}' for v in tmp_vcfs)
                    + f' --OUTPUT {output_vcf}'
                ),
                input_files_or_dirs=tmp_vcfs,
                output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
            )
            self.remove_files_and_dirs(
                *chain.from_iterable(
                    [o.path for o in t] for t in input_targets
                )
            )


class CNNScoreVariants(SagvcTask):
    input_vcf_path = luigi.Parameter()
    input_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    interval_list_path = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    gatk = luigi.Parameter(default='gatk')
    python = luigi.Parameter(default='python')
    add_cnnscorevariants_args = luigi.ListParameter(default=list())
    save_memory = luigi.BoolParameter(default=False)
    message = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.vcf.gz' + s)
            for s in ['', '.tbi']
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        input_vcf = Path(self.input_vcf_path).resolve()
        input_cram = Path(self.input_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        interval_list = Path(self.interval_list_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        output_vcf = output_files[0]
        self.setup_shell(
            run_id='.'.join(output_vcf.name.split('.')[:-2]),
            commands=[self.gatk, self.python], cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} CNNScoreVariants'
                + f' --input {input_cram}'
                + f' --variant {input_vcf}'
                + f' --reference {fa}'
                + f' --intervals {interval_list}'
                + ''.join(
                    f' {a}' for a in [
                        *self.add_cnnscorevariants_args,
                        *(
                            ['--disable-bam-index-caching', 'true']
                            if self.save_memory else list()
                        )
                    ]
                )
                + f' --output {output_vcf}'
            ),
            input_files_or_dirs=[
                input_vcf, fa, input_cram, interval_list
            ],
            output_files_or_dirs=output_files
        )


@requires(ScoreVariantsWithCnn)
class FilterVariantTranches(SagvcTask):
    resource_vcf_paths = luigi.ListParameter()
    gatk = luigi.Parameter(default='gatk')
    add_filtervarianttranches_args = luigi.ListParameter(
        default=[
            '--info-key', 'CNN_2D', '--snp-tranche', '99.9',
            '--snp-tranche', '99.95', '--indel-tranche', '99.0',
            '--indel-tranche', '99.4', '--invalidate-previous-filters'
        ]
    )
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_path_prefix = re.sub(r'\.vcf\.gz$', '', self.input()[0][0].path)
        return [
            luigi.LocalTarget(f'{output_path_prefix}.filtered.vcf.gz{s}')
            for s in ['', '.tbi']
        ]

    def run(self):
        input_vcf = Path(self.input()[0].path)
        run_id = '.'.join(input_vcf.name.split('.')[:-3])
        self.print_log(f'Apply tranche filtering:\t{run_id}')
        resource_vcfs = [Path(p) for p in self.resource_vcf_paths]
        output_files = [Path(o.path) for o in self.output()]
        output_vcf = output_files[0]
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} FilterVariantTranches'
                + f' --variant {input_vcf}'
                + ''.join(f' --resource {p}' for p in resource_vcfs)
                + ''.join(
                    f' {a}' for a in [
                        *self.add_filtervarianttranches_args,
                        *(
                            ['--disable-bam-index-caching', 'true']
                            if self.save_memory else list()
                        )
                    ]
                )
                + f' --output {output_vcf}'
            ),
            input_files_or_dirs=[input_vcf, *resource_vcfs],
            output_files_or_dirs=output_files
        )


if __name__ == '__main__':
    luigi.run()
