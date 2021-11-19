#!/usr/bin/env python

import os
from itertools import product
from math import floor
from pathlib import Path

import luigi

from .core import SagvcTask


class CallSomaticStructualVariantsWithManta(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    python2 = luigi.Parameter(default='python2')
    configmantapy_path = luigi.Parameter(default='./configManta.py')
    exome = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        tn_stem = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{tn_stem}.manta.{v}SV.vcf.gz{s}')
            ) for v, s in product(['somatic', 'diploid'], ['', '.tbi'])
        ]

    def run(self):
        output_links = [Path(o.path) for o in self.output()]
        run_id = '.'.join(output_links[0].name.split('.')[:-4])
        self.print_log(f'Call somatic SVs with Manta:\t{run_id}')
        tumor_cram = Path(self.tumor_cram_path).resolve()
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        bed = (Path(self.bed_path).resolve() if self.bed_path else None)
        config_script = Path(self.configmantapy_path).resolve()
        dest_dir = output_links[0].parent
        run_dir = dest_dir.joinpath(run_id)
        run_script = run_dir.joinpath('runWorkflow.py')
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(
                [
                    'somaticSV', 'diploidSV', 'candidateSV',
                    'candidateSmallIndels'
                ],
                ['', '.tbi']
            )
        ]
        pythonpath = '{0}:{1}'.format(
            config_script.parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = max(floor(self.memory_mb / 1024), 4)
        self.setup_shell(
            run_id=run_id, commands=[self.python2, config_script],
            cwd=run_dir, **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {config_script}'
                + f' --tumorBam={tumor_cram}'
                + f' --normalBam={normal_cram}'
                + f' --referenceFasta={fa}'
                + f' --runDir={run_dir}'
                + (f' --callRegions={bed}' if bed else '')
                + (' --exome' if self.exome else '')
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa,
                *([bed] if self.bed_path else list())
            ],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[
                run_script, tumor_cram, normal_cram, fa,
                *([bed] if self.bed_path else list())
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for o in output_links:
            f = run_dir.joinpath('results/variants').joinpath(
                o.name.split('.manta.')[-1]
            ).relative_to(dest_dir)
            self.run_shell(args=f'ln -s {f} {o}', output_files_or_dirs=o)


class CallGermlineStructualVariantsWithManta(SagvcTask):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    python2 = luigi.Parameter(default='python2')
    configmantapy_path = luigi.Parameter(default='./configManta.py')
    exome = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        cram_stem = Path(self.normal_cram_path).stem
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{cram_stem}.manta.diploidSV.vcf.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_links = [Path(o.path) for o in self.output()]
        run_id = '.'.join(output_links[0].name.split('.')[:-4])
        self.print_log(f'Call germline SVs with Manta:\t{run_id}')
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        bed = (Path(self.bed_path).resolve() if self.bed_path else None)
        config_script = Path(self.configmantapy_path).resolve()
        dest_dir = output_links[0].parent
        run_dir = dest_dir.joinpath(run_id)
        run_script = run_dir.joinpath('runWorkflow.py')
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(
                ['diploidSV', 'candidateSV', 'candidateSmallIndels'],
                ['', '.tbi']
            )
        ]
        pythonpath = '{0}:{1}'.format(
            config_script.parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = max(floor(self.memory_mb / 1024), 4)
        self.setup_shell(
            run_id=run_id, commands=[self.python2, config_script],
            cwd=run_dir, **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {config_script}'
                + f' --bam={normal_cram}'
                + f' --referenceFasta={fa}'
                + f' --runDir={run_dir}'
                + (f' --callRegions={bed}' if bed else '')
                + (' --exome' if self.exome else '')
            ),
            input_files_or_dirs=[
                normal_cram, fa, *([bed] if self.bed_path else list())
            ],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[
                run_script, normal_cram, fa,
                *([bed] if self.bed_path else list())
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for o in output_links:
            f = run_dir.joinpath('results/variants').joinpath(
                o.name.split('.manta.')[-1]
            ).relative_to(dest_dir)
            self.run_shell(args=f'ln -s {f} {o}', output_files_or_dirs=o)


if __name__ == '__main__':
    luigi.run()
