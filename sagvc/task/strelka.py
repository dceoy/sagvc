#!/usr/bin/env python

import os
from itertools import product
from math import floor
from pathlib import Path

import luigi

from .core import SagvcTask
from .manta import CallSomaticStructualVariantsWithManta


class CallSomaticVariantsWithStrelka(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    manta_indel_vcf_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    python2 = luigi.Parameter(default='python2')
    bcftools = luigi.Parameter(default='bcftools')
    configurestrelkasomaticworkflowpy_path = luigi.Parameter(
        default='./configureStrelkaSomaticWorkflow.py'
    )
    configmantapy_path = luigi.Parameter(default='./configManta.py')
    exome = luigi.BoolParameter(default=False)
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def require(self):
        if self.manta_indel_vcf_path:
            return super().requires()
        else:
            return CallSomaticStructualVariantsWithManta(
                tumor_cram_path=self.tumor_cram_path,
                normal_cram_path=self.normal_cram_path, fa_path=self.fa_path,
                bed_path=self.bed_path, dest_dir_path=self.dest_dir_path,
                python2=self.python2,
                configmantapy_path=self.configmantapy_path, exome=self.exome,
                n_cpu=self.n_cpu, memory_mb=self.memory_mb,
                sh_config=self.sh_config
            )

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        tn_stem = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{tn_stem}.strelka.somatic.vcf.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_id = '.'.join(output_vcf.name.split('.')[:-4])
        self.print_log(f'Call somatic variants with Strelka:\t{run_id}')
        tumor_cram = Path(self.tumor_cram_path).resolve()
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        bed = (Path(self.bed_path).resolve() if self.bed_path else None)
        manta_indel_vcf = Path(
            self.manta_indel_vcf_path
            or Path(self.input()[0].path).parent.joinpath(
                f'{run_id}/results/variants/candidateSmallIndels.vcf.gz'
            )
        ).resolve()
        config_script = Path(
            self.configurestrelkasomaticworkflowpy_path
        ).resolve()
        dest_dir = output_vcf.parent
        run_dir = dest_dir.joinpath(run_id)
        run_script = run_dir.joinpath('runWorkflow.py')
        result_files = [
            run_dir.joinpath(f'results/variants/somatic.{v}.vcf.gz{s}')
            for v, s in product(['snvs', 'indels'], ['', '.tbi'])
        ]
        pythonpath = '{0}:{1}'.format(
            Path(config_script).parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = floor(self.memory_mb / 1024)
        self.setup_shell(
            run_id=run_id,
            commands=[self.python2, config_script, self.bcftools],
            cwd=run_dir, **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {config_script}'
                + f' --tumorBam={tumor_cram}'
                + f' --normalBam={normal_cram}'
                + f' --referenceFasta={fa}'
                + f' --indelCandidates={manta_indel_vcf}'
                + f' --callRegions={bed}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.exome else '')
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa, manta_indel_vcf, bed
            ],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[
                run_script, tumor_cram, normal_cram, fa, manta_indel_vcf, bed
            ],
            output_files_or_dirs=[*result_files, run_dir]
        )
        self.bcftools_concat(
            input_vcf_paths=[
                o for o in result_files if o.name.endswith('.vcf.gz')
            ],
            output_vcf_path=output_vcf, bcftools=self.bcftools,
            n_cpu=self.n_cpu, memory_mb=self.memory_mb, index_vcf=True,
            remove_input=False
        )


class CallGermlineVariantsWithStrelka(SagvcTask):
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    python2 = luigi.Parameter(default='python2')
    bcftools = luigi.Parameter(default='bcftools')
    configurestrelkagermlineworkflowpy_path = luigi.Parameter(
        default='./configureStrelkaGermlineWorkflow.py'
    )
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
                dest_dir.joinpath(f'{cram_stem}.strelka.germline.vcf.gz{s}')
            ) for s in ['', '.tbi']
        ]

    def run(self):
        output_links = [Path(o.path) for o in self.output()]
        run_id = '.'.join(output_links[0].name.split('.')[:-4])
        self.print_log(f'Call germline variants with Strelka:\t{run_id}')
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        bed = (Path(self.bed_path).resolve() if self.bed_path else None)
        config_script = Path(self.cf['configureStrelkaGermlineWorkflow.py'])
        config_script = Path(
            self.configurestrelkagermlineworkflowpy_path
        ).resolve()
        dest_dir = output_links[0].parent
        run_dir = dest_dir.joinpath(run_id)
        run_script = run_dir.joinpath('runWorkflow.py')
        result_files = [
            run_dir.joinpath(f'results/variants/{v}.vcf.gz{s}')
            for v, s in product(['variants', 'genome'], ['', '.tbi'])
        ]
        pythonpath = '{0}:{1}'.format(
            Path(config_script).parent.parent.joinpath('lib/python'),
            (os.getenv('PYTHONPATH') or '')
        )
        memory_gb = floor(self.memory_mb / 1024)
        self.setup_shell(
            run_id=run_id, commands=[self.python2, config_script], cwd=run_dir,
            **self.sh_config, env={'PYTHONPATH': pythonpath}
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {config_script}'
                + f' --bam={normal_cram}'
                + f' --referenceFasta={fa}'
                + f' --callRegions={bed}'
                + f' --runDir={run_dir}'
                + (' --exome' if self.exome else '')
            ),
            input_files_or_dirs=[normal_cram, fa, bed],
            output_files_or_dirs=[run_script, run_dir]
        )
        self.run_shell(
            args=(
                f'set -e && {self.python2} {run_script} --mode=local'
                + f' --jobs={self.n_cpu} --memGb={memory_gb}'
            ),
            input_files_or_dirs=[run_script, normal_cram, fa, bed],
            output_files_or_dirs=[*result_files, run_dir]
        )
        for o in output_links:
            f = run_dir.joinpath('results/variants').joinpath(
                'variants.' + o.name.split('.strelka.germline.')[-1]
            ).relative_to(dest_dir)
            self.run_shell(args=f'ln -s {f} {o}', output_files_or_dirs=o)


if __name__ == '__main__':
    luigi.run()
