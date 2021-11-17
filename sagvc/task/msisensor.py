#!/usr/bin/env python

from pathlib import Path

import luigi
from ftarc.task.samtools import SamtoolsView

from .core import SagvcTask


class ScanMicrosatellites(SagvcTask):
    fa_path = luigi.Parameter()
    dest_dir_path = luigi.Parameter(default='.')
    msisensor = luigi.Parameter(default='msisensor')
    add_scan_args = luigi.ListParameter(default=list())
    sh_config = luigi.DictParameter(default=dict())
    priority = 60

    def output(self):
        return luigi.LocalTarget(
            Path(self.dest_dir_path).resolve().joinpath(
                Path(self.fa_path).stem + '.microsatellites.tsv'
            )
        )

    def run(self):
        run_id = Path(self.fa_path).stem
        self.print_log(f'Scan microsatellites:\t{run_id}')
        fa = Path(self.fa_path).resolve()
        output_tsv = Path(self.output().path)
        self.setup_shell(
            run_id=run_id, commands=self.msisensor, cwd=output_tsv.parent,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.msisensor} scan'
                + f' -d {fa}'
                + ''.join(f' {a}' for a in self.add_scan_args)
                + f' -o {output_tsv}'
            ),
            input_files_or_dirs=fa, output_files_or_dirs=output_tsv
        )


class ScoreMsiWithMsisensor(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    microsatellites_tsv_path = luigi.Parameter(default='')
    bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    msisensor = luigi.Parameter(default='msisensor')
    samtools = luigi.Parameter(default='samtools')
    add_msi_args = luigi.ListParameter(default=list())
    n_cpu = luigi.IntParameter(default=1)
    sh_config = luigi.DictParameter(default=dict())
    priority = 20

    def requires(self):
        if self.microsatellites_tsv_path:
            return super().requires()
        else:
            return ScanMicrosatellites(
                fa_path=self.fa_path,
                dest_dir_path=str(
                    Path(self.dest_dir_path).joinpath(
                        Path(self.fa_path).stem + '.wgs'
                    )
                ),
                msisensor=self.msisensor, sh_config=self.sh_config
            )

    def output(self):
        dest_dir = Path(self.dest_dir_path).resolve()
        tn_stem = self.create_matched_id(
            self.tumor_cram_path, self.normal_cram_path
        )
        return [
            luigi.LocalTarget(
                dest_dir.joinpath(f'{tn_stem}.msisensor.tsv{s}')
            ) for s in ['', '_dis', '_germline', '_somatic']
        ]

    def run(self):
        output_files = [Path(i.path) for i in self.output()]
        run_dir = output_files[0].parent
        crams = [
            Path(p).resolve()
            for p in [self.tumor_cram_path, self.normal_cram_path]
        ]
        input_targets = yield [
            SamtoolsView(
                input_sam_path=str(c),
                output_sam_path=str(run_dir.joinpath(f'{c.stem}.bam')),
                fa_path=self.fa_path, samtools=self.samtools, n_cpu=self.n_cpu,
                remove_input=False, index_sam=True, sh_config=self.sh_config
            ) for c in crams
        ]
        bams = [Path(i[0].path) for i in input_targets]
        run_id = Path(output_files[0].stem).stem
        self.print_log(f'Score MSI with MSIsensor:\t{run_id}')
        ms_tsv = Path(
            self.microsatellites_tsv_path or self.input().path
        ).resource()
        bed = (Path(self.bed_path).resolve() if self.bed_path else None)
        self.setup_shell(
            run_id=run_id, commands=self.msisensor, cwd=run_dir,
            **self.sh_config
        )
        self.run_shell(
            args=(
                f'set -e && {self.msisensor} msi'
                + f' -d {ms_tsv}'
                + (' -e {bed}' if bed else '')
                + ''.join(f' {a}' for a in self.add_msi_args)
                + f' -o {output_files[0]}'
                + f' -t {bams[0]} -n {bams[1]}'
            ),
            input_files_or_dirs=[*bams, ms_tsv, *([bed] if bed else list())],
            output_files_or_dirs=[*output_files, run_dir]
        )
        for c, t in zip(crams, input_targets):
            if str(c) != t[0].path:
                self.remove_files_and_dirs(*[i.path for i in t])


if __name__ == '__main__':
    luigi.run()
