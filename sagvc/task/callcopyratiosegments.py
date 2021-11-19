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
    add_preprocessintervals_args = luigi.ListParameter(
        default=['--interval-merging-rule', 'OVERLAPPING_ONLY']
    )
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.interval_list_path or self.fa_path).stem + (
                    ('.excl.' + Path(self.cnv_blacklist_path).stem)
                    if self.cnv_blacklist_path else ''
                ) + '.w{}s.preproc.interval_list'.format(
                    'x' if self.exome else 'g'
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
        add_preprocessintervals_args = [
            *self.add_preprocessintervals_args,
            *(
                list() if '--bin-length' in self.add_preprocessintervals_args
                else ['--bin-length', (0 if self.exome else 1000)]
            ),
            *(
                list() if '--padding' in self.add_preprocessintervals_args
                else ['--padding', (250 if self.exome else 0)]
            )
        ]
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
                + ''.join(f' {a}' for a in add_preprocessintervals_args)
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
    add_collectalleliccounts_args = luigi.ListParameter(default=list())
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.cram_path).stem + '.allelic_counts.tsv'
            )
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
                + f' --input {cram}'
                + f' --reference {fa}'
                + f' --intervals {snp_interval_list}'
                + ''.join(
                    f' {a}' for a in [
                        *self.add_collectalleliccounts_args,
                        *(
                            ['--disable-bam-index-caching', 'true']
                            if self.save_memory else list()
                        )
                    ]
                )
                + f' --output {allelic_counts_tsv}'
            ),
            input_files_or_dirs=[cram, snp_interval_list, fa],
            output_files_or_dirs=allelic_counts_tsv
        )


class CollectReadCounts(SagvcTask):
    cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    preproc_interval_list_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    exome = luigi.BoolParameter(default=False)
    add_preprocessintervals_args = luigi.ListParameter(
        default=['--interval-merging-rule', 'OVERLAPPING_ONLY']
    )
    add_collectreadcounts_args = luigi.ListParameter(
        default=['--interval-merging-rule', 'OVERLAPPING_ONLY']
    )
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def requires(self):
        if self.preproc_interval_list_path:
            return super().requires()
        else:
            return PreprocessIntervals(
                fa_path=self.fa_path,
                cnv_blacklist_path=self.cnv_blacklist_path,
                interval_list_path=self.interval_list_path,
                dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                exome=self.exome,
                add_preprocessintervals_args=self.add_preprocessintervals_args,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.cram_path).stem + '.counts.hdf5'
            )
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
                + ''.join(
                    f' {a}' for a in [
                        *self.add_collectreadcounts_args,
                        *(
                            ['--disable-bam-index-caching', 'true']
                            if self.save_memory else list()
                        )
                    ]
                )
                + f' --output {counts_hdf5}'
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
    add_denoisereadcounts_args = luigi.ListParameter(default=list())
    add_plotdenoisedcopyratios_args = luigi.ListParameter(default=list())
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
        fa_dict = Path(self.fa_dict_path).resolve()
        dest_dir = denoised_cr_tsv.parent
        self.setup_shell(
            run_id=run_id, commands=[self.gatk, self.r], cwd=dest_dir,
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
                + ''.join(f' {a}' for a in self.add_denoisereadcounts_args)
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
                    + ''.join(
                        f' {a}' for a in self.add_plotdenoisedcopyratios_args
                    )
                    + f' --output {dest_dir}'
                    + f' --output-prefix {run_id}'

                ),
                input_files_or_dirs=[
                    standardized_cr_tsv, denoised_cr_tsv, fa_dict
                ],
                output_files_or_dirs=dest_dir.joinpath(
                    f'{run_id}.denoised.png'
                )
            )


class ModelSegments(SagvcTask):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    fa_dict_path = luigi.Parameter()
    snp_interval_list_path = luigi.Parameter()
    tumor_cram_path = luigi.Parameter(default='')
    preproc_interval_list_path = luigi.Parameter(default='')
    cnv_blacklist_path = luigi.Parameter(default='')
    interval_list_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    r = luigi.Parameter(default='R')
    exome = luigi.BoolParameter(default=False)
    create_plots = luigi.BoolParameter(default=True)
    add_preprocessintervals_args = luigi.ListParameter(
        default=['--interval-merging-rule', 'OVERLAPPING_ONLY']
    )
    add_collectreadcounts_args = luigi.ListParameter(
        default=['--interval-merging-rule', 'OVERLAPPING_ONLY']
    )
    add_denoisereadcounts_args = luigi.ListParameter(default=list())
    add_plotdenoisedcopyratios_args = luigi.ListParameter(default=list())
    add_modelsegments_args = luigi.ListParameter(default=list())
    add_plotmodeledsegments_args = luigi.ListParameter(default=list())
    save_memory = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def requires(self):
        return [
            DenoiseReadCounts(
                cram_path=(self.tumor_cram_path or self.normal_cram_path),
                fa_path=self.fa_path, fa_dict_path=self.fa_dict_path,
                preproc_interval_list_path=self.preproc_interval_list_path,
                cnv_blacklist_path=self.cnv_blacklist_path,
                interval_list_path=self.interval_list_path,
                dest_dir_path=self.dest_dir_path, gatk=self.gatk, r=self.r,
                exome=self.exome, create_plots=self.create_plots,
                add_preprocessintervals_args=self.add_preprocessintervals_args,
                add_collectreadcounts_args=self.add_collectreadcounts_args,
                add_denoisereadcounts_args=self.add_denoisereadcounts_args,
                add_plotdenoisedcopyratios_args=(
                    self.add_plotdenoisedcopyratios_args
                ),
                save_memory=self.save_memory, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ),
            *[
                CollectAllelicCounts(
                    cram_path=p, fa_path=self.fa_path,
                    snp_interval_list_path=self.snp_interval_list_path,
                    dest_dir_path=self.dest_dir_path, gatk=self.gatk,
                    save_memory=self.save_memory, n_cpu=self.n_cpu,
                    memory_mb=self.memory_mb, sh_config=self.sh_config
                ) for p in [self.tumor_cram_path, self.normal_cram_path] if p
            ]
        ]

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        output_stem = (
            self.create_matched_id(
                self.tumor_cram_path, self.normal_cram_path
            ) if self.tumor_cram_path else Path(self.normal_cram_path).stem
        )
        return [
            luigi.LocalTarget(dest_dir.joinpath(f'{output_stem}.{s}'))
            for s in ['cr.seg', 'hets.tsv', 'modelFinal.seg']
        ]

    def run(self):
        output_files = [Path(o.path) for o in self.output()]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Produce denoised copy ratios:\t{run_id}')
        denoised_cr_tsv = Path(self.input()[0][0].path)
        normal_allelic_counts_tsv = Path(self.input()[-1].path)
        dest_dir = output_files[0].parent
        if self.tumor_cram_path:
            case_allelic_counts_tsv = Path(self.input()[1].path)
            input_files = [
                denoised_cr_tsv, case_allelic_counts_tsv,
                normal_allelic_counts_tsv
            ]
            add_modelsegments_args = [
                *self.add_modelsegments_args,
                '--allelic-counts', case_allelic_counts_tsv,
                '--normal-allelic-counts', normal_allelic_counts_tsv,
                *(
                    list() if '--minimum-total-allele-count-case'
                    in self.add_modelsegments_args
                    else ['--minimum-total-allele-count-case', 0]
                ),
                *(
                    list() if '--minimum-total-allele-count-normal'
                    in self.add_modelsegments_args
                    else ['--minimum-total-allele-count-normal', 30]
                )
            ]
        else:
            input_files = [denoised_cr_tsv, normal_allelic_counts_tsv]
            add_modelsegments_args = [
                *self.add_modelsegments_args,
                '--allelic-counts', normal_allelic_counts_tsv,
                *(
                    list() if '--minimum-total-allele-count-case'
                    in self.add_modelsegments_args
                    else ['--minimum-total-allele-count-case', 30]
                )
            ]
        self.setup_shell(
            run_id=run_id, commands=[self.gatk, self.r], cwd=dest_dir,
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
                + ''.join(f' {a}' for a in add_modelsegments_args)
                + f' --output {dest_dir}'
                + f' --output-prefix {run_id}'
            ),
            input_files_or_dirs=input_files,
            output_files_or_dirs=output_files
        )
        if self.create_plots:
            fa_dict = Path(self.fa_dict_path).resolve()
            het_allelic_counts_tsv = output_files[1]
            modeled_segments = output_files[2]
            self.run_shell(
                args=(
                    f'set -e && {self.gatk} PlotModeledSegments'
                    + f' --denoised-copy-ratios {denoised_cr_tsv}'
                    + f' --allelic-counts {het_allelic_counts_tsv}'
                    + f' --segments {modeled_segments}'
                    + f' --sequence-dictionary {fa_dict}'
                    + ''.join(
                        f' {a}' for a in self.add_plotmodeledsegments_args
                    )
                    + f' --output {dest_dir}'
                    + f' --output-prefix {run_id}'
                ),
                input_files_or_dirs=[
                    denoised_cr_tsv, het_allelic_counts_tsv, modeled_segments,
                    fa_dict
                ],
                output_files_or_dirs=dest_dir.joinpath(f'{run_id}.modeled.png')
            )


@requires(ModelSegments)
class CallCopyRatioSegments(SagvcTask):
    dest_dir_path = luigi.Parameter(default='.')
    gatk = luigi.Parameter(default='gatk')
    add_callcopyratiosegments_args = luigi.ListParameter(default=list())
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
                + f' --input {cr_seg}'
                + ''.join(f' {a}' for a in self.add_callcopyratiosegments_args)
                + f' --output {output_seg}'
            ),
            input_files_or_dirs=cr_seg, output_files_or_dirs=output_seg
        )


@requires(CallCopyRatioSegments)
class CallCopyRatioSegmentsTumor(luigi.WrapperTask):
    priority = 30

    def output(self):
        return self.input()


@requires(CallCopyRatioSegments)
class CallCopyRatioSegmentsNormal(luigi.WrapperTask):
    tumor_cram_path = ''
    priority = 30

    def output(self):
        return self.input()


@requires(CallCopyRatioSegmentsTumor, CallCopyRatioSegmentsNormal)
class CallCopyRatioSegmentsMatched(luigi.WrapperTask):
    priority = 30

    def output(self):
        return self.input()


if __name__ == '__main__':
    luigi.run()
