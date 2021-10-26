#!/usr/bin/env python

from pathlib import Path

import luigi

from .core import SagvcTask


class CallSomaticStructualVariantsWithDelly(SagvcTask):
    tumor_cram_path = luigi.Parameter()
    normal_cram_path = luigi.Parameter()
    fa_path = luigi.Parameter()
    tumor_sample_name = luigi.Parameter()
    normal_sample_name = luigi.Parameter()
    excl_bed_path = luigi.Parameter(default='')
    dest_dir_path = luigi.Parameter(default='.')
    delly = luigi.Parameter(default='delly')
    bcftools = luigi.Parameter(default='bcftools')
    n_cpu = luigi.IntParameter(default=1)
    memory_mb = luigi.FloatParameter(default=4096)
    sh_config = luigi.DictParameter(default=dict())
    priority = 30

    def output(self):
        run_dir = Path(self.dest_dir_path).resolve().joinpath(
            self.create_matched_id(self.tumor_cram_path, self.normal_cram_path)
        )
        return [
            luigi.LocalTarget(
                run_dir.joinpath(run_dir.name + f'.delly.{s}')
            ) for s in [
                'filtered.vcf.gz', 'filtered.vcf.gz.tbi', 'filtered.bcf',
                'filtered.bcf.csi', 'bcf', 'bcf.csi'
            ]
        ]

    def run(self):
        output_vcf = Path(self.output()[0].path)
        run_dir = output_vcf.parent
        run_id = run_dir.name
        self.print_log(f'Call somatic SVs with Delly:\t{run_id}')
        input_crams = [Path(self.tumor_cram_path), Path(self.normal_cram_path)]
        fa = Path(self.fa_path)
        excl_bed = (Path(self.excl_bed_path) if self.excl_bed_path else None)
        raw_bcf = Path(self.output()[4].path)
        filtered_bcf = Path(self.output()[2].path)
        samples_tsv = run_dir.joinpath(f'{raw_bcf.stem}.tsv')
        self.setup_shell(
            run_id=run_id, commands=[self.delly, self.bcftools], cwd=run_dir,
            **self.sh_config, env={'OMP_NUM_THREADS': str(self.n_cpu)}
        )
        self.run_shell(
            args=(
                f'set -e && {self.delly} call'
                + f' --genome {fa}'
                + (f' --exclude {excl_bed}' if excl_bed else '')
                + f' --outfile {raw_bcf}'
                + ''.join(f' {p}' for p in input_crams)
            ),
            input_files_or_dirs=[
                *input_crams, fa,
                *([excl_bed] if self.excl_bed_path else list())
            ],
            output_files_or_dirs=[raw_bcf, f'{raw_bcf}.csi', run_dir]
        )
        self.run_shell(
            args=(
                'set -e && echo -ne "'
                + f'{self.tumor_sample_name}\\ttumor\\n'
                + f'{self.normal_sample_name}\\tcontrol\\n'
                + f'" | tee {samples_tsv}'
            ),
            output_files_or_dirs=samples_tsv
        )
        self.run_shell(
            args=(
                f'set -e && {self.delly} filter --filter somatic'
                + f' --samples {samples_tsv}'
                + f' --outfile {filtered_bcf}'
                + f' {raw_bcf}'
            ),
            input_files_or_dirs=[raw_bcf, samples_tsv],
            output_files_or_dirs=[filtered_bcf, f'{filtered_bcf}.csi']
        )
        self.bcftools_sort(
            input_vcf_path=filtered_bcf, output_vcf_path=output_vcf,
            bcftools=self.bcftools, n_cpu=self.n_cpu, memory_mb=self.memory_mb,
            index_vcf=True, remove_input=False
        )


if __name__ == '__main__':
    luigi.run()
