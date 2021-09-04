#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.resource import FetchReferenceFasta
from luigi.util import requires

from .core import SagvcTask
from .samtools import SamtoolsView


@requires(FetchReferenceFasta)
class ScanMicrosatellites(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    msisensorpro = luigi.Parameter(default='msisensor-pro')
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        fa = Path(self.input()[0].path)
        return luigi.LocalTarget(
            fa.parent.joinpath(f'{fa.stem}.microsatellites.tsv')
        )

    def run(self):
        fa = Path(self.input()[0].path)
        run_id = fa.stem
        self.print_log(f'Scan microsatellites:\t{run_id}')
        output_tsv = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.msisensorpro, cwd=self.dest_dir_path,
            **self.sh_config
        )
        self.run_shell(
            args=f'set -e && {self.msisensorpro} scan -d {fa} -o {output_tsv}',
            input_files_or_dirs=fa, output_files_or_dirs=output_tsv
        )


@requires(FetchReferenceFasta)
class ScoreMsiWithMsisensorPro(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    samtools = luigi.Parameter(default='samtools')
    microsatellites_tsv_path = luigi.Parameter()
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    msisensorpro = luigi.Parameter(default='msisensor-pro')
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        output_tsv_name = '{}.msisensorpro.tsv'.format(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return [
            luigi.LocalTarget(dest_dir.joinpath(output_tsv_name + s))
            for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        bam_targets = yield [
            SamtoolsView(
                input_sam_path=p,
                output_sam_path=str(dest_dir.joinpath(Path(p).stem + '.bam')),
                fa_path=self.fa_path, samtools=self.samtools, n_cpu=self.n_cpu,
                remove_input=False, index_sam=True, sh_config=self.sh_config
            ) for p in [self.tumor_cram_path, self.normal_cram_path]
        ]
        output_files = [Path(i.path) for i in self.input()]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Score MSI with MSIsensor-pro:\t{run_id}')
        bams = [Path(i[0].path) for i in bam_targets]
        microsatellites_tsv = Path(self.microsatellites_tsv_path).resource()
        bed = (Path(self.bed_path).resource() if self.bed_path else None)
        self.setup_shell(
            run_id=run_id, commands=self.msisensorpro, cwd=dest_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.msisensorpro} msi'
                + f' -t {bams[0]} -n {bams[1]}'
                + f' -d {microsatellites_tsv}'
                + (' -e {bed}' if bed else '')
                + f' -o {output_files[0]}'
            ),
            input_files_or_dirs=[
                *bams, microsatellites_tsv, *([bed] if bed else list())
            ],
            output_files_or_dirs=output_files
        )
        self.remove_files_and_dirs(
            *[i.path for i in bam_targets[0]],
            *[i.path for i in bam_targets[1]]
        )


if __name__ == '__main__':
    luigi.run()
