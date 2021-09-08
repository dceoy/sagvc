#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta

from .core import SagvcTask
from .resource import CreateRegionListBed


class CreateCnvAcccessBed(SagvcTask):
    fa_path = luigi.Parameter()
    cnv_blacklist_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    pigz = luigi.Parameter(default='pigz')
    pbzip2 = luigi.Parameter(default='pbzip2')
    samtools = luigi.Parameter(default='samtools')
    gatk = luigi.Parameter(default='gatk')
    bgzip = luigi.Parameter(default='bgzip')
    tabix = luigi.Parameter(default='tabix')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def requires(self):
        return [
            FetchReferenceFasta(
                fa_path=self.fa_path, pigz=self.pigz, pbzip2=self.pbzip2,
                samtools=self.samtools, gatk=self.gatk, n_cpu=self.n_cpu,
                memory_mb=self.memory_mb, sh_config=self.sh_config
            ),
            CreateRegionListBed(
                region_list_path=self.cnv_blacklist_path,
                dest_dir_path=self.dest_dir_path, bgzip=self.bgzip,
                tabix=self.tabix, n_cpu=self.n_cpu, sh_config=self.sh_config
            )
        ]

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.fa_path).stem + '.not_in.'
                + Path(self.cnv_blacklist_path).stem + '.access.bed'
            )
        )

    def run(self):
        fa = Path(self.input()[0][0].path)
        run_id = fa.stem
        self.print_log(
            f'Calculate accessible coordinates for CNV calling:\t{run_id}'
        )
        cnv_blacklist_bed = Path(self.input()[1][-1].path)
        output_bed = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.cnvkitpy, cwd=output_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} access'
                + f' --exclude={cnv_blacklist_bed}'
                + f' --output={output_bed} {fa}'
            ),
            input_files_or_dirs=[fa, cnv_blacklist_bed],
            output_files_or_dirs=output_bed
        )


class CallSomaticCnvWithCnvkit(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    refflat_txt_path = luigi.Parameter()
    access_bed_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    samtools = luigi.Parameter(default='samtools')
    rscript = luigi.Parameter(default='Rscript')
    seq_method = luigi.Parameter(default='wgs')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        tumor_stem = Path(self.tumor_cram_path).stem
        normal_stem = Path(self.normal_cram_path).stem
        return [
            luigi.LocalTarget(dest_dir.joinpath(n)) for n in (
                [
                    (tumor_stem + s) for s in [
                        '.call.seg', '.call.cns', '.cns', '.bintest.cns',
                        '.cnr', '.targetcoverage.cnn',
                        '.antitargetcoverage.cnn', '-diagram.pdf',
                        '-scatter.png'
                    ]
                ] + [
                    (normal_stem + s) for s in [
                        '.targetcoverage.cnn', '.antitargetcoverage.cnn',
                        '.reference.cnn'
                    ]
                ]
            )
        ]

    def run(self):
        run_id = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        self.print_log(f'Score MSI with MSIsensor-pro:\t{run_id}')
        tumor_cram = Path(self.tumor_cram_path).resolve()
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        access_bed = Path(self.access_bed).resolve()
        refflat_txt = Path(self.refflat_txt_path).resolve()
        output_files = [Path(o.path) for o in self.output()]
        output_call_cns = output_files[0]
        dest_dir = output_call_cns.parent
        output_ref_cnn = dest_dir.joinpath(f'{normal_cram.stem}.reference.cnn')
        output_call_seg = dest_dir.joinpath(f'{output_call_cns.stem}.seg')
        self.setup_shell(
            run_id=run_id,
            commands=[self.cnvkitpy, self.samtools, self.rscript],
            cwd=dest_dir, **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} batch'
                + f' --seq-method={self.seq_method}'
                + f' --fasta={fa}'
                + f' --access={access_bed}'
                + f' --annotate={refflat_txt}'
                + f' --processes={self.n_cpu}'
                + ' --drop-low-coverage --diagram --scatter'
                + f' --output-dir={dest_dir}'
                + f' --output-reference={output_ref_cnn}'
                + f' --normal={normal_cram}'
                + f' {tumor_cram}'
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa, access_bed, refflat_txt
            ],
            output_files_or_dirs=output_files[1:]
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} export seg'
                + f' --output={output_call_seg}'
                + f' {output_call_cns}'
            ),
            input_files_or_dirs=output_call_cns,
            output_files_or_dirs=output_call_seg
        )


if __name__ == '__main__':
    luigi.run()
