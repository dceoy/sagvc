#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.samtools import SamtoolsView

from .core import SagvcTask


class ScanMicrosatellites(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    msisensor = luigi.Parameter(default='msisensor')
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.fa_path).stem + '.microsatellites.tsv'
            )
        )

    def run(self):
        fa = Path(self.input()[0].path)
        run_id = fa.stem
        self.print_log(f'Scan microsatellites:\t{run_id}')
        output_tsv = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.msisensor, cwd=output_tsv.parent,
            **self.sh_config
        )
        self.run_shell(
            args=f'set -e && {self.msisensor} scan -d {fa} -o {output_tsv}',
            input_files_or_dirs=fa, output_files_or_dirs=output_tsv
        )


class ScoreMsiWithMsisensor(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    microsatellites_tsv_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    msisensor = luigi.Parameter(default='msisensor')
    samtools = luigi.Parameter(default='samtools')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        output_tsv_name = '{}.msisensor.tsv'.format(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return [
            luigi.LocalTarget(dest_dir.joinpath(output_tsv_name + s))
            for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        crams = [
            Path(p) for p in [self.tumor_cram_path, self.normal_cram_path]
        ]
        bams = [
            (c if c.suffix == '.bam' else dest_dir.joinpath(f'{c.stem}.bam'))
            for c in crams
        ]
        tmp_target = yield [
            SamtoolsView(
                input_sam_path=str(c), output_sam_path=str(b),
                fa_path=self.fa_path, samtools=self.samtools, n_cpu=self.n_cpu,
                remove_input=False, index_sam=True, sh_config=self.sh_config
            ) for c, b in zip(crams, bams) if c != b
        ]
        output_files = [Path(i.path) for i in self.output()]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Score MSI with MSIsensor:\t{run_id}')
        microsatellites_tsv = Path(self.microsatellites_tsv_path).resource()
        bed = (Path(self.bed_path).resource() if self.bed_path else None)
        self.setup_shell(
            run_id=run_id, commands=self.msisensor, cwd=dest_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.msisensor} msi'
                + f' -t {bams[0]} -n {bams[1]}'
                + f' -d {microsatellites_tsv}'
                + (' -e {bed}' if bed else '')
                + f' -o {output_files[0]}'
            ),
            input_files_or_dirs=[
                *bams, microsatellites_tsv,
                *([bed] if bed else list())
            ],
            output_files_or_dirs=output_files
        )
        for t in tmp_target:
            self.remove_files_and_dirs(*[i.path for i in t])


if __name__ == '__main__':
    luigi.run()
