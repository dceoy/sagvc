#!/usr/bin/env python

from pathlib import Path

import luigi

from .core import SagvcTask


class CreateCnvAcccessBed(SagvcTask):
    fa_path = luigi.Parameter()
    excl_bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    add_access_args = luigi.ListParameter(default=list())
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 10

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.fa_path).stem + (
                    ('.excl.' + Path(self.excl_bed_path).stem)
                    if self.excl_bed_path else ''
                ) + '.access.bed'
            )
        )

    def run(self):
        run_id = Path(self.fa_path).stem
        self.print_log(
            f'Calculate accessible coordinates for CNV calling:\t{run_id}'
        )
        fa = Path(self.fa_path).resolve()
        excl_bed = (
            Path(self.excl_bed_path).resolve() if self.excl_bed_path else None
        )
        output_bed = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.cnvkitpy, cwd=output_bed.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} access'
                + (f' --exclude={excl_bed}' if excl_bed else '')
                + ''.join(f' {a}' for a in self.add_access_args)
                + f' --output={output_bed} {fa}'
            ),
            input_files_or_dirs=[fa, *([excl_bed] if excl_bed else list())],
            output_files_or_dirs=output_bed
        )


class CallSomaticCnvWithCnvkit(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    refflat_txt_path = luigi.Parameter(default='')
    access_bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    cnvkitpy = luigi.Parameter(default='cnvkit.py')
    samtools = luigi.Parameter(default='samtools')
    rscript = luigi.Parameter(default='Rscript')
    seq_method = luigi.Parameter(default='wgs')
    diagram = luigi.BoolParameter(default=True)
    scatter = luigi.BoolParameter(default=True)
    add_batch_args = luigi.ListParameter(
        default=['--drop-low-coverage', ' --short-names']
    )
    add_export_seg_args = luigi.ListParameter(default=list())
    add_diagram_args = luigi.ListParameter(default=list())
    add_scatter_args = luigi.ListParameter(default=list())
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        tumor_stem = Path(self.tumor_cram_path).stem
        normal_stem = Path(self.normal_cram_path).stem
        access_stem = Path(self.access_bed_path or self.fa_path).stem
        return [
            luigi.LocalTarget(run_dir.joinpath(n)) for n in (
                [
                    (tumor_stem + s) for s in (
                        [
                            '.call.seg', '.seg', '.call.cns', '.cns',
                            '.bintest.cns', '.cnr', '.targetcoverage.cnn',
                            '.antitargetcoverage.cnn'
                        ] + (['-diagram.pdf'] if self.diagram else list())
                        + (['-scatter.pdf'] if self.scatter else list())
                    )
                ] + [
                    (normal_stem + s) for s in [
                        '.targetcoverage.cnn', '.antitargetcoverage.cnn',
                        '.reference.cnn'
                    ]
                ] + [
                    (access_stem + s) for s in [
                        '.target.bed', '.antitarget.bed'
                    ]
                ]
            )
        ]

    def run(self):
        run_id = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        self.print_log(f'Call somatic CNVs with CNVkit:\t{run_id}')
        tumor_cram = Path(self.tumor_cram_path).resolve()
        normal_cram = Path(self.normal_cram_path).resolve()
        fa = Path(self.fa_path).resolve()
        access_bed = (
            Path(self.access_bed_path).resolve()
            if self.access_bed_path else None
        )
        refflat_txt = (
            Path(self.refflat_txt_path).resolve()
            if self.refflat_txt_path else None
        )
        output_files = [Path(o.path) for o in self.output()]
        run_dir = output_files[0].parent
        output_ref_cnn = run_dir.joinpath(f'{normal_cram.stem}.reference.cnn')
        output_call_cns = output_files[2]
        output_cns = output_files[3]
        self.setup_shell(
            run_id=run_id,
            commands=[self.cnvkitpy, self.samtools, self.rscript], cwd=run_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.cnvkitpy} batch'
                + f' --seq-method={self.seq_method}'
                + f' --fasta={fa}'
                + (f' --access={access_bed}' if access_bed else '')
                + (f' --annotate={refflat_txt}' if refflat_txt else '')
                + f' --output-dir={run_dir}'
                + f' --output-reference={output_ref_cnn}'
                + f' --processes={self.n_cpu}'
                + ''.join(f' {a}' for a in self.add_batch_args)
                + f' --normal={normal_cram}'
                + f' {tumor_cram}'
            ),
            input_files_or_dirs=[
                tumor_cram, normal_cram, fa,
                *[f for f in [access_bed, refflat_txt] if f]
            ],
            output_files_or_dirs=[
                *[f for f in output_files[2:] if f.suffix != '.pdf'], run_dir
            ]
        )
        for o in [output_call_cns, output_cns]:
            output_seg = run_dir.joinpath(f'{o.stem}.seg')
            self.run_shell(
                args=(
                    f'set -e && {self.cnvkitpy} export seg'
                    + ''.join(f' {a}' for a in self.add_export_seg_args)
                    + f' --output={output_seg} {o}'
                ),
                input_files_or_dirs=o, output_files_or_dirs=output_seg
            )
        for c in ['diagram', 'scatter']:
            if getattr(self, c):
                graph_pdf = run_dir.joinpath(f'{output_cns.stem}.{c}.pdf')
                self.run_shell(
                    args=(
                        f'set -e && {self.cnvkitpy} {c}'
                        + ''.join(
                            f' {a}' for a in getattr(self, f'add_{c}_args')
                        )
                        + f' --output={graph_pdf} {output_cns}'
                    ),
                    input_files_or_dirs=output_cns,
                    output_files_or_dirs=graph_pdf
                )


if __name__ == '__main__':
    luigi.run()
