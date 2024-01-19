from os import makedirs
from os.path import basename, join
from typing import List, Optional
from .parse_vcf import ParseVcf
from .template import Settings, Processor
from .variant_picking import VariantPicking
from .variant_filtering import FlagVariants, RemoveVariants


def variant_filtering(
        input_vcf: str,
        output_vcf: str,
        variant_flagging_criteria: str,
        variant_removal_flags: str,
        only_pass: bool,
        workdir: str):

    makedirs(workdir, exist_ok=True)

    settings = Settings(
        workdir=workdir,
        outdir='.',
        threads=1,
        debug=False,
        mock=False)

    VariantFiltering(settings).main(
        input_vcf=input_vcf,
        output_vcf=output_vcf,
        variant_flagging_criteria=variant_flagging_criteria,
        variant_removal_flags=variant_removal_flags,
        only_pass=only_pass)


class VariantFiltering(Processor):

    input_vcf: str
    output_vcf: str
    variant_flagging_criteria: str
    variant_removal_flags: List[str]
    only_pass: bool

    vcf: str

    def main(
            self,
            input_vcf: str,
            output_vcf: str,
            variant_flagging_criteria: str,
            variant_removal_flags: str,
            only_pass: bool):

        self.input_vcf = input_vcf
        self.output_vcf = output_vcf
        self.variant_flagging_criteria = variant_flagging_criteria
        self.variant_removal_flags = [] if variant_removal_flags.lower() == 'none' else variant_removal_flags.split(',')
        self.only_pass = only_pass

        self.vcf = self.input_vcf
        self.flag_variants()
        self.remove_variants()
        self.call(f'mv {self.vcf} {output_vcf}')

    def flag_variants(self):
        self.vcf = FlagVariants(self.settings).main(
            vcf=self.input_vcf,
            variant_flagging_criteria=self.variant_flagging_criteria)

    def remove_variants(self):
        self.vcf = RemoveVariants(self.settings).main(
            vcf=self.vcf,
            flags=self.variant_removal_flags,
            only_pass=self.only_pass)


def variant_picking(
        ref_fa: str,
        output_vcf: str,
        mutect2: str,
        haplotype_caller: str,
        muse: str,
        lofreq: str,
        varscan: str,
        vardict: str,
        somatic_sniper: str,
        min_snv_callers: int,
        min_indel_callers: int,
        workdir: str):

    makedirs(workdir, exist_ok=True)

    settings = Settings(
        workdir=workdir,
        outdir='.',
        threads=1,
        debug=False,
        mock=False)

    Picking(settings).main(
        ref_fa=ref_fa,
        output_vcf=output_vcf,
        mutect2=None if mutect2.lower() == 'none' else mutect2,
        haplotype_caller=None if haplotype_caller.lower() == 'none' else haplotype_caller,
        muse=None if muse.lower() == 'none' else muse,
        lofreq=None if lofreq.lower() == 'none' else lofreq,
        varscan=None if varscan.lower() == 'none' else varscan,
        vardict=None if vardict.lower() == 'none' else vardict,
        somatic_sniper=None if somatic_sniper.lower() == 'none' else somatic_sniper,
        min_snv_callers=min_snv_callers,
        min_indel_callers=min_indel_callers
    )


class Picking(Processor):

    ref_fa: str
    output_vcf: str
    mutect2: Optional[str]
    haplotype_caller: Optional[str]
    muse: Optional[str]
    lofreq: Optional[str]
    varscan: Optional[str]
    vardict: Optional[str]
    somatic_sniper: Optional[str]
    min_snv_callers: int
    min_indel_callers: int

    vcfs: List[str]

    def main(
            self,
            ref_fa: str,
            output_vcf: str,
            mutect2: Optional[str],
            haplotype_caller: Optional[str],
            muse: Optional[str],
            lofreq: Optional[str],
            varscan: Optional[str],
            vardict: Optional[str],
            somatic_sniper: Optional[str],
            min_snv_callers: int,
            min_indel_callers: int):

        self.ref_fa = ref_fa
        self.output_vcf = output_vcf
        self.mutect2 = mutect2
        self.haplotype_caller = haplotype_caller
        self.muse = muse
        self.lofreq = lofreq
        self.varscan = varscan
        self.vardict = vardict
        self.somatic_sniper = somatic_sniper
        self.min_snv_callers = min_snv_callers
        self.min_indel_callers = min_indel_callers

        self.copy_vcfs()
        self.pick_variants()

    def copy_vcfs(self):
        self.vcfs = []
        for src, dst in [
            (self.mutect2, f'{self.workdir}/mutect2.vcf'),
            (self.haplotype_caller, f'{self.workdir}/haplotype-caller.vcf'),
            (self.muse, f'{self.workdir}/muse.vcf'),
            (self.lofreq, f'{self.workdir}/lofreq.vcf'),
            (self.varscan, f'{self.workdir}/varscan.vcf'),
            (self.vardict, f'{self.workdir}/vardict.vcf'),
            (self.somatic_sniper, f'{self.workdir}/somatic-sniper.vcf'),
        ]:
            if src is not None:
                if src.endswith('.gz'):
                    dst += '.gz'
                self.call(f'cp {src} {dst}')
                self.vcfs.append(dst)

    def pick_variants(self):
        vcf = VariantPicking(self.settings).main(
            ref_fa=self.ref_fa,
            vcfs=self.vcfs,
            min_snv_callers=self.min_snv_callers,
            min_indel_callers=self.min_indel_callers)

        self.call(f'mv {vcf} {self.output_vcf}')


def vcf2csv(
        input_vcf: str,
        output_csv: str,
        workdir: str):

    makedirs(workdir, exist_ok=True)

    settings = Settings(
        workdir=workdir,
        outdir='.',
        threads=1,
        debug=False,
        mock=False)

    Vcf2Csv(settings).main(
        input_vcf=input_vcf,
        output_csv=output_csv)


class Vcf2Csv(Processor):

    input_vcf: str
    output_csv: str

    def main(
            self,
            input_vcf: str,
            output_csv: str):

        self.input_vcf = input_vcf
        self.output_csv = output_csv

        ParseVcf(self.settings).main(
            vcf=self.input_vcf,
            dstdir=self.workdir)

        csv = join(self.workdir, basename(self.input_vcf)[:-len('.vcf')] + '.csv')
        self.call(f'mv {csv} {self.output_csv}')
