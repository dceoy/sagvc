#!/usr/bin/env python

import re
from pathlib import Path

import luigi
from luigi.util import requires

from .core import SagvcTask


class PreprocessIntervals(SagvcTask):
    fa_path = luigi.Parameter()
    cnv_blacklist_path = luigi.Parameter(default='')
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    exome = luigi.BoolParameter(default=False)
    param_dict = luigi.DictParameter(default=dict())
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                '{0}.preprocessed.{1}.interval_list'.format(
                    '.'.join([
                        Path(p).stem
                        for p in [self.fa_path, self.interval_list_path] if p
                    ]),
                    ('wxs' if self.exome else 'wgs')
                )
            )
        )

    def run(self):
        output_interval_list = Path(self.output().path)
        run_id = output_interval_list.stem
        self.print_log(f'Prepare bins for coverage collection:\t{run_id}')
        fa = Path(self.fa_path).resolve()
        interval_list = (
            Path(self.interval_list_path).resolve()
            if self.interval_list_path else None
        )
        cnv_blacklist = (
            Path(self.cnv_blacklist_path).resolve()
            if self.cnv_blacklist_path else None
        )
        param_dict = (
            self.param_dict or (
                {'bin-length': 0, 'padding': 250} if self.exome
                else {'bin-length': 1000, 'padding': 0}
            )
        )
        self.setup_shell(
            run_id=run_id, commands=self.gatk,
            cwd=output_interval_list.parent, **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} PreprocessIntervals'
                + f' --reference {fa}'
                + (f' --intervals {interval_list}' if interval_list else '')
                + (
                    f' --exclude-intervals {cnv_blacklist}'
                    if cnv_blacklist else ''
                )
                + ''.join(f' --{k} {v}' for k, v in param_dict.items())
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {output_interval_list}'
            ),
            input_files_or_dirs=[
                fa, *[p for p in [interval_list, cnv_blacklist] if p]
            ],
            output_files_or_dirs=output_interval_list
        )


class CollectAllelicCounts(SagvcTask):
    cram_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.cram_path).stem
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.allelic_counts.tsv')
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Collects allele counts:\t{run_id}')
        allelic_counts_tsv = Path(self.output().path)
        cram = Path(self.cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        snp_interval_list = Path(self.snp_interval_list_path).resolve()
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=allelic_counts_tsv.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} CollectAllelicCounts'
                + f' --intervals {snp_interval_list}'
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --output {allelic_counts_tsv}'
                + ' --disable-bam-index-caching '
                + str(self.save_memory).lower()
            ),
            input_files_or_dirs=[cram, snp_interval_list, fa],
            output_files_or_dirs=allelic_counts_tsv
        )


class CollectAllelicCountsTumor(luigi.WrapperTask):
    tumor_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def requires(self):
        return CollectAllelicCounts(
            cram_path=self.tumor_cram_path, fa_path=self.fa_path,
            snp_interval_list_path=self.snp_interval_list_path,
            dest_dir_path=self.dest_dir_path, gatk=self.gatk,
            save_memory=self.save_memory, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class CollectAllelicCountsNormal(luigi.WrapperTask):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def requires(self):
        return CollectAllelicCounts(
            cram_path=self.normal_cram_path, fa_path=self.fa_path,
            snp_interval_list_path=self.snp_interval_list_path,
            dest_dir_path=self.dest_dir_path, gatk=self.gatk,
            save_memory=self.save_memory, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, sh_config=self.sh_config
        )

    def output(self):
        return self.input()


class CollectReadCounts(SagvcTask):
    cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    preproc_interval_list_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.cram_path).stem
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.counts.hdf5')
        )

    def run(self):
        run_id = Path(self.cram_path).stem
        self.print_log(f'Collects read counts:\t{run_id}')
        cram = Path(self.cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        preproc_interval_list = Path(self.preproc_interval_list_path).resolve()
        counts_hdf5 = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=counts_hdf5.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} CollectReadCounts'
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --intervals {preproc_interval_list}'
                + ' --format HDF5'
                + ' --interval-merging-rule OVERLAPPING_ONLY'
                + f' --output {counts_hdf5}'
                + ' --disable-bam-index-caching '
                + str(self.save_memory).lower()
            ),
            input_files_or_dirs=[cram, fa, preproc_interval_list],
            output_files_or_dirs=counts_hdf5
        )


@requires(CollectReadCounts)
class DenoiseReadCounts(SagvcTask):
    fa_dict_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    r = luigi.Parameter(default='R')
    create_plots = luigi.BoolParameter(default=True)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        output_path_prefix = re.sub(r'\.counts\.hdf5$', '', self.input().path)
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}')
            for s in ['denoised_cr.tsv', 'standardized_cr.tsv']
        ]

    def run(self):
        counts_hdf5 = Path(self.input().path)
        run_id = Path(counts_hdf5.stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        denoised_cr_tsv = Path(self.output()[0].path)
        standardized_cr_tsv = Path(self.output()[1].path)
        fa_dict = Path(self.fa_dict_path)
        run_dir = denoised_cr_tsv.parent
        self.setup_shell(
            run_id=run_id, commands=[self.gatk, self.r], cwd=run_dir,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} DenoiseReadCounts'
                + f' --input {counts_hdf5}'
                + f' --standardized-copy-ratios {standardized_cr_tsv}'
                + f' --denoised-copy-ratios {denoised_cr_tsv}'
            ),
            input_files_or_dirs=counts_hdf5,
            output_files_or_dirs=[standardized_cr_tsv, denoised_cr_tsv]
        )
        if self.create_plots:
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} PlotDenoisedCopyRatios'
                    + f' --standardized-copy-ratios {standardized_cr_tsv}'
                    + f' --denoised-copy-ratios {denoised_cr_tsv}'
                    + f' --sequence-dictionary {fa_dict}'
                    + f' --output {run_dir}'
                    + f' --output-prefix {run_id}'

                ),
                input_files_or_dirs=[
                    standardized_cr_tsv, denoised_cr_tsv, fa_dict
                ],
                output_files_or_dirs=run_dir.joinpath(f'{run_id}.denoised.png')
            )


@requires(DenoiseReadCounts)
class ModelSegments(SagvcTask):
    normal_allelic_counts_tsv_path = luigi.Parameter()
    fa_dict_path = luigi.Parameter()
    case_allelic_counts_tsv_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    r = luigi.Parameter(default='R')
    dest_dir_path = luigi.Parameter(default='.')
    min_total_allele_counts = luigi.DictParameter(
        default={'case': 0, 'normal': 30}
    )
    create_plots = luigi.BoolParameter(default=True)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        output_path_prefix = str(
            Path(self.dest_dir_path).joinpath(
                Path(
                    self.create_matched_id(
                        self.case_allelic_counts_tsv_path,
                        self.normal_allelic_counts_tsv_path
                    ) if self.case_allelic_counts_tsv_path else
                    Path(self.normal_allelic_counts_tsv_path).stem
                ).stem
            )
        )
        return [
            luigi.LocalTarget(f'{output_path_prefix}.{s}')
            for s in ['cr.seg', 'hets.tsv', 'modelFinal.seg']
        ]

    def run(self):
        output_files = [Path(o.path) for o in self.output()]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        denoised_cr_tsv = Path(self.input()[0].path)
        normal_allelic_counts_tsv = Path(
            self.normal_allelic_counts_tsv_path
        ).resolve()
        run_dir = output_files[0].parent
        if self.case_allelic_counts_tsv_path:
            case_allelic_counts_tsv = Path(
                self.case_allelic_counts_tsv_path
            ).resolve()
            allelic_count_args = (
                f' --allelic-counts {case_allelic_counts_tsv}'
                + f' --normal-allelic-counts {normal_allelic_counts_tsv}'
                + ''.join(
                    f' --minimum-total-allele-count-{k} {v}'
                    for k, v in self.min_total_allele_counts.items()
                )
            )
            input_files = [
                denoised_cr_tsv, case_allelic_counts_tsv,
                normal_allelic_counts_tsv
            ]
        else:
            allelic_count_args = (
                f' --allelic-counts {normal_allelic_counts_tsv}'
                + ' --minimum-total-allele-count-case {}'.format(
                    self.min_total_allele_counts['normal']
                )
            )
            input_files = [denoised_cr_tsv, normal_allelic_counts_tsv]
        self.setup_shell(
            run_id=run_id, commands=[self.gatk, self.r], cwd=run_dir,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }

        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} ModelSegments'
                + f' --denoised-copy-ratios {denoised_cr_tsv}'
                + allelic_count_args
                + f' --output-prefix {run_id}'
                + f' --output {run_dir}'
            ),
            input_files_or_dirs=input_files,
            output_files_or_dirs=[*output_files, run_dir]
        )
        if self.create_plots:
            plots_dir = run_dir.joinpath(f'{run_id}.plots')
            fa_dict = Path(self.fa_dict_path)
            het_allelic_counts_tsv = output_files[1]
            modeled_segments = output_files[2]
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} PlotModeledSegments'
                    + f' --denoised-copy-ratios {denoised_cr_tsv}'
                    + f' --allelic-counts {het_allelic_counts_tsv}'
                    + f' --segments {modeled_segments}'
                    + f' --sequence-dictionary {fa_dict}'
                    + f' --output {plots_dir}'
                    + f' --output-prefix {run_id}'
                ),
                input_files_or_dirs=[
                    denoised_cr_tsv, het_allelic_counts_tsv, modeled_segments,
                    fa_dict
                ],
                output_files_or_dirs=[
                    plots_dir.joinpath(f'{run_id}.modeled.png'), plots_dir
                ]
            )


@requires(ModelSegments)
class CallCopyRatioSegments(SagvcTask):
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            re.sub(r'\.seg$', '.called.seg', self.input()[0].path)
        )

    def run(self):
        cr_seg = Path(self.input()[0].path)
        run_id = cr_seg.stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        output_seg = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.gatk, cwd=output_seg.parent,
            **self.sh_config,
            env={
                'JAVA_TOOL_OPTIONS': self.generate_gatk_java_options(
                    n_cpu=self.n_cpu, memory_mb=self.memory_mb
                )
            }
        )
        self.run_shell(
            args=(
                f'set -e && {self.gatk} CallCopyRatioSegments'
                + f' --input {cr_seg} --output {output_seg}'
            ),
            input_files_or_dirs=cr_seg, output_files_or_dirs=output_seg
        )


class CallCopyRatioSegmentsTumor(luigi.Task):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    preproc_interval_list_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    r = luigi.Parameter(default='R')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.cr.called.seg')
        )

    def run(self):
        input_targets = [
            *[
                CollectAllelicCounts(
                    cram_path=p, fa_path=self.fa_path,
                    snp_interval_list_path=self.snp_interval_list_path,
                    dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                    save_memory=self.save_memory, n_cpu=self.n_cpu,
                    memory_mb=self.memory_mb, sh_config=self.sh_config
                ) for p in [self.tumor_cram_path, self.normal_cram_path]
            ],
            *(
                list() if self.preproc_interval_list_path else [
                    PreprocessIntervals(
                        fa_path=self.fa_path,
                        cnv_blacklist_path=self.cnv_blacklist_path,
                        interval_list_path=self.interval_list_path,
                        dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                        n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                        sh_config=self.sh_config
                    )
                ]
            )
        ]
        yield CallCopyRatioSegments(
            cram_path=self.tumor_cram_path, fa_path=self.fa_path,
            fa_dict_path=str(
                Path(self.fa_path).parent.joinpath(
                    Path(self.fa_path).stem + '.dict'
                )
            ),
            case_allelic_counts_tsv_path=input_targets[0].path,
            normal_allelic_counts_tsv_path=input_targets[1].path,
            preproc_interval_list_path=(
                self.preproc_interval_list_path or self.input()[2].path
            ),
            dest_dir_path=self.dest_dir_path, gatk=self.gatk, r=self.r,
            save_memory=self.save_memory, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, sh_config=self.sh_config
        )


class CallCopyRatioSegmentsNormal(luigi.Task):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    preproc_interval_list_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    r = luigi.Parameter(default='R')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            Path(self.normal_cram_path).stem
        )
        return luigi.LocalTarget(
            run_dir.joinpath(f'{run_dir.name}.cr.called.seg')
        )

    def run(self):
        input_targets = [
            CollectAllelicCounts(
                cram_path=self.normal_cram_path, fa_path=self.fa_path,
                snp_interval_list_path=self.snp_interval_list_path,
                dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ),
            *(
                list() if self.preproc_interval_list_path else [
                    PreprocessIntervals(
                        fa_path=self.fa_path,
                        cnv_blacklist_path=self.cnv_blacklist_path,
                        interval_list_path=self.interval_list_path,
                        dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                        n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                        sh_config=self.sh_config
                    )
                ]
            )
        ]
        yield CallCopyRatioSegments(
            cram_path=self.normal_cram_path, fa_path=self.fa_path,
            fa_dict_path=str(
                Path(self.fa_path).parent.joinpath(
                    Path(self.fa_path).stem + '.dict'
                )
            ),
            case_allelic_counts_tsv_path='',
            normal_allelic_counts_tsv_path=input_targets[0].path,
            preproc_interval_list_path=(
                self.preproc_interval_list_path or self.input()[1].path
            ),
            dest_dir_path=self.dest_dir_path, gatk=self.gatk, r=self.r,
            save_memory=self.save_memory, n_cpu=self.n_cpu,
            memory_mb=self.memory_mb, sh_config=self.sh_config
        )


@requires(CallCopyRatioSegmentsTumor, CallCopyRatioSegmentsNormal)
class CallCopyRatioSegmentsMatched(luigi.WrapperTask):
    priority = 30

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
