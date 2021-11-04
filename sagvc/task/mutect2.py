#!/usr/bin/env python

import re
from itertools import chain
from pathlib import Path

import luigi
from luigi.util import requires

from .core import SagvcTask
from .resource import SplitEvaluationIntervals


class GetPileupSummaries(SagvcTask):
    cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    interval_list_path = luigi.Parameter()
    common_biallelic_vcf_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.cram_path).stem
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.pileup.table')
        )

    def run(self):
        cram = Path(self.cram_path).resolve()
        run_id = cram.stem
        self.print_log(f'Get pileup summary:\t{run_id}')
        output_pileup_table = Path(self.output().path)
        fa = Path(self.fa_path).resolve()
        interval_list = Path(self.interval_list_path).resolve()
        common_biallelic_vcf = Path(
            self.common_biallelic_vcf_path
        ).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_pileup_table.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} GetPileupSummaries'
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --variant {common_biallelic_vcf}'
                + f' --intervals {interval_list}'
                + f' --output {output_pileup_table}'
                + ' --disable-bam-index-caching '
                + str(self.save_memory).lower()
            ),
            input_files_or_dirs=[
                cram, fa, interval_list, common_biallelic_vcf
            ],
            output_files_or_dirs=output_pileup_table
        )


class CalculateContamination(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    common_biallelic_vcf_path = luigi.Parameter()
    interval_list_path = luigi.Parameter(default='')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return [
            luigi.LocalTarget(run_dir.joinpath(f'{run_dir.name}.{s}.table'))
            for s in ['contamination', 'segment']
        ]

    def run(self):
        output_contamination_table = Path(self.output()[0].path)
        run_dir = output_contamination_table.parent
        input_targets = yield [
            GetPileupSummaries(
                cram_path=p, fa_path=self.fa_path,
                interval_list_path=self.interval_list_path,
                common_biallelic_vcf_path=self.common_biallelic_vcf_path,
                dest_dir_path=str(run_dir), gatk=self.gatk,
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ) for p in [self.tumor_cram_path, self.normal_cram_path]
        ]
        run_id = '.'.join(output_contamination_table.name.split('.')[:-2])
        self.print_log(f'Calculate cross-sample contamination:\t{run_id}')
        pileup_tables = [Path(i.path) for i in input_targets]
        output_segment_table = Path(self.output()[1].path)
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
                f'set -e && {self.gatk} CalculateContamination'
                + f' --input {pileup_tables[0]}'
                + f' --matched-normal {pileup_tables[1]}'
                + f' --output {output_contamination_table}'
                + f' --tumor-segmentation {output_segment_table}'
            ),
            input_files_or_dirs=pileup_tables,
            output_files_or_dirs=[
                output_contamination_table, output_segment_table
            ]
        )


@requires(SplitEvaluationIntervals)
class CallVariantsWithMutect2(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    tumor_sample_name = luigi.Parameter()
    normal_sample_name = luigi.Parameter()
    germline_resource_vcf_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    samtools = luigi.Parameter(default='samtools')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return [
            luigi.LocalTarget(run_dir.joinpath(f'{run_dir.name}.mutect2.{s}'))
            for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'cram', 'cram.crai',
                'read-orientation-model.tar.gz'
            ]
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
            Mutect2(
                input_cram_paths=[self.tumor_cram_path, self.normal_cram_path],
                fa_path=str(fa),
                germline_resource_vcf_path=self.germline_resource_vcf_path,
                interval_list_path=str(o),
                tumor_sample_name=self.tumor_sample_name,
                normal_sample_name=self.normal_sample_name,
                output_path_prefix=s, gatk=self.gatk,
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ) for o, s in zip(interval_lists, tmp_prefixes)
        ]
        run_id = '.'.join(output_vcf.name.split('.')[:-3])
        self.print_log(f'Call somatic variants with Mutect2:\t{run_id}')
        output_stats = Path(self.output()[2].path)
        output_cram = Path(self.output()[3].path)
        ob_priors = Path(self.output()[5].path)
        f1r2s = [f'{s}.f1r2.tar.gz' for s in tmp_prefixes]
        self.setup_shell(
            run_id=run_id, commands=[self.gatk, self.samtools],
            cwd=output_vcf.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} LearnReadOrientationModel'
                + ''.join(f' --input {f}' for f in f1r2s)
                + f' --output {ob_priors}'
            ),
            input_files_or_dirs=f1r2s, output_files_or_dirs=ob_priors
        )
        if skip_interval_list_split:
            tmp_bam = f'{tmp_prefixes[0]}.bam'
            self.samtools_view(
                input_sam_path=tmp_bam, fa_path=fa,
                output_sam_path=output_cram, samtools=self.samtools,
                n_cpu=self.n_cpu, index_sam=True, remove_input=True
            )
        else:
            tmp_vcfs = [Path(f'{s}.vcf.gz') for s in tmp_prefixes]
            self.picard_mergevcfs(
                input_vcf_paths=tmp_vcfs, output_vcf_path=output_vcf,
                picard=self.gatk, remove_input=False
            )
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} MergeVcfs'
                    + ''.join(f' --INPUT {v}' for v in tmp_vcfs)
                    + f' --REFERENCE_SEQUENCE {fa}'
                    + f' --OUTPUT {output_vcf}'
                ),
                input_files_or_dirs=tmp_vcfs,
                output_files_or_dirs=[output_vcf, f'{output_vcf}.tbi']
            )
            tmp_statses = [Path(f'{s}.vcf.gz.stats') for s in tmp_prefixes]
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} MergeMutectStats'
                    + ''.join(f' --stats {s}' for s in tmp_statses)
                    + f' --output {output_stats}'
                ),
                input_files_or_dirs=tmp_statses,
                output_files_or_dirs=output_stats
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


class Mutect2(SagvcTask):
    input_cram_paths = luigi.ListParameter()
    fa_path = luigi.Parameter()
    interval_list_path = luigi.Parameter()
    tumor_sample_name = luigi.Parameter()
    normal_sample_name = luigi.Parameter()
    output_path_prefix = luigi.Parameter()
    germline_resource_vcf_path = luigi.Parameter(default='')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    message = luigi.Parameter(default='')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        return [
            luigi.LocalTarget(f'{self.output_path_prefix}.{s}') for s in [
                'vcf.gz', 'vcf.gz.tbi', 'vcf.gz.stats', 'bam', 'f1r2.tar.gz'
            ]
        ]

    def run(self):
        if self.message:
            self.print_log(self.message)
        input_crams = [Path(p).resolve() for p in self.input_cram_paths]
        fa = Path(self.fa_path).resolve()
        germline_resource_vcf = (
            Path(self.germline_resource_vcf_path).resolve()
            if self.germline_resource_vcf_path else None
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
                f'set -e && {self.gatk} Mutect2'
                + ''.join(f' --input {c}' for c in input_crams)
                + f' --reference {fa}'
                + f' --intervals {interval_list}'
                + (
                    f' --germline-resource {germline_resource_vcf}'
                    if germline_resource_vcf else ''
                )
                + f' --output {output_vcf}'
                + f' --bam-output {output_files[3]}'
                + f' --f1r2-tar-gz {output_files[4]}'
                + f' --tumor-sample {self.tumor_sample_name}'
                + f' --normal-sample {self.normal_sample_name}'
                + f' --native-pair-hmm-threads {self.n_cpu}'
                + ' --max-mnp-distance 0'
                + ' --create-output-bam-index false'
                + ' --disable-bam-index-caching '
                + str(self.save_memory).lower()
            ),
            input_files_or_dirs=[
                *input_crams, fa, interval_list,
                *([germline_resource_vcf] if germline_resource_vcf else list())
            ],
            output_files_or_dirs=[*output_files, run_dir]
        )


@requires(CallVariantsWithMutect2, CalculateContamination)
class FilterMutectCalls(SagvcTask):
    fa_path = luigi.Parameter()
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 50

    def output(self):
        output_vcf_path = re.sub(
            r'\.vcf\.gz$', '.filtered.vcf.gz', self.input()[0][0].path
        )
        return [
            luigi.LocalTarget(output_vcf_path + s)
            for s in ['', '.tbi', '.stats']
        ]

    def run(self):
        input_vcf = Path(self.input()[0][0].path)
        run_id = '.'.join(input_vcf.name.split('.')[:-3])
        self.print_log(f'Filter somatic variants called by Mutect2:\t{run_id}')
        input_stats = Path(self.input()[0][2].path)
        ob_priors = Path(self.input()[0][5].path)
        fa = Path(self.fa_path).resolve()
        contamination_table = Path(self.input()[1][0].path)
        segment_table = Path(self.input()[1][1].path)
        output_files = [Path(o.path) for o in self.output()]
        output_vcf = output_files[0]
        output_stats = output_files[2]
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
                f'set -e && {self.gatk} FilterMutectCalls'
                + f' --reference {fa}'
                + f' --variant {input_vcf}'
                + f' --stats {input_stats}'
                + f' --contamination-table {contamination_table}'
                + f' --tumor-segmentation {segment_table}'
                + f' --orientation-bias-artifact-priors {ob_priors}'
                + f' --output {output_vcf}'
                + f' --filtering-stats {output_stats}'
            ),
            input_files_or_dirs=[
                input_vcf, fa, input_stats, ob_priors, contamination_table,
                segment_table,
            ],
            output_files_or_dirs=output_files
        )


if __name__ == '__main__':
    luigi.run()