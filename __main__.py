import argparse
from typing import List, Dict
from src import variant_filtering, variant_picking, vcf2csv, remove_umi


__VERSION__ = '1.2.1-beta'


PROG = 'python omic'
DESCRIPTION = f'CLI tools for variant calling pipeline (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'


VARIANT_FILTERING = 'variant-filtering'
VARIANT_PICKING = 'variant-picking'
VCF2CSV = 'vcf2csv'
REMOVE_UMI = 'remove-umi'


INPUT_VCF_ARG = {
    'keys': ['-i', '--input-vcf'],
    'properties': {
        'type': str,
        'required': True,
        'help': 'path to the input vcf(.gz) file',
    }
}
OUTPUT_VCF_ARG = {
    'keys': ['-o', '--output-vcf'],
    'properties': {
        'type': str,
        'required': True,
        'help': 'path to the output vcf(.gz) file',
    }
}
WORKDIR_ARG = {
    'keys': ['-w', '--workdir'],
    'properties': {
        'type': str,
        'required': False,
        'default': './temp',
        'help': 'path to the temporary working directory (default: %(default)s)',
    }
}
HELP_ARG = {
    'keys': ['-h', '--help'],
    'properties': {
        'action': 'help',
        'help': 'show this help message',
    }
}
VERSION_ARG = {
    'keys': ['-v', '--version'],
    'properties': {
        'action': 'version',
        'version': __VERSION__,
        'help': 'show version',
    }
}


MODE_TO_GROUP_TO_ARGS = {
    VARIANT_FILTERING:
        {
            'Required':
                [
                    INPUT_VCF_ARG,
                    OUTPUT_VCF_ARG
                ],
            'Optional':
                [
                    {
                        'keys': ['--variant-flagging-criteria'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'comma-separated flagging criteria, e.g. "low_depth: DP<20, mid_qual: 20<=MQ<=40, high_af: AF>0.02" (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--variant-removal-flags'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'comma-separated flags for variant removal, e.g. "panel_of_normals,map_qual" (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--only-pass'],
                        'properties': {
                            'action': 'store_true',
                            'help': 'only keep the variants with PASS in FILTER column',
                        }
                    },
                    WORKDIR_ARG,
                    HELP_ARG,
                    VERSION_ARG,
                ],
        },
    VARIANT_PICKING:
        {
            'Required':
                [
                    {
                        'keys': ['-r', '--ref-fa'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the reference genome fasta file',
                        }
                    },
                    OUTPUT_VCF_ARG,
                ],
            'Optional':
                [
                    {
                        'keys': ['--mutect2'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the mutect2 vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--haplotype-caller'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the haplotype-caller vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--muse'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the muse vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--varscan'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the varscan vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--lofreq'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the lofreq vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--vardict'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the vardict vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--somatic-sniper'],
                        'properties': {
                            'type': str,
                            'required': False,
                            'default': 'None',
                            'help': 'path to the somatic-sniper vcf(.gz) file (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--min-snv-callers'],
                        'properties': {
                            'type': int,
                            'required': False,
                            'default': 1,
                            'help': 'min number of variant callers for an SNV to be picked (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['--min-indel-callers'],
                        'properties': {
                            'type': int,
                            'required': False,
                            'default': 1,
                            'help': 'min number of variant callers for an indel to be picked (default: %(default)s)',
                        }
                    },
                    WORKDIR_ARG,
                    HELP_ARG,
                    VERSION_ARG,
                ],
        },
    VCF2CSV:
        {
            'Required':
                [
                    INPUT_VCF_ARG,
                    {
                        'keys': ['-o', '--output-csv'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the output csv file',
                        }
                    },
                ],
            'Optional':
                [
                    WORKDIR_ARG,
                    HELP_ARG,
                    VERSION_ARG,
                ],
        },
    REMOVE_UMI:
        {
            'Required':
                [
                    {
                        'keys': ['-1', '--input-fq1'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the input read 1 fastq(.gz) file',
                        }
                    },
                    {
                        'keys': ['-2', '--input-fq2'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the input read 2 fastq(.gz) file',
                        }
                    },
                    {
                        'keys': ['-3', '--output-fq1'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the output read 1 fastq(.gz) file',
                        }
                    },
                    {
                        'keys': ['-4', '--output-fq2'],
                        'properties': {
                            'type': str,
                            'required': True,
                            'help': 'path to the output read 2 fastq(.gz) file',
                        }
                    },
                ],
            'Optional':
                [
                    {
                        'keys': ['-l', '--umi-length'],
                        'properties': {
                            'type': int,
                            'required': False,
                            'default': 0,
                            'help': 'UMI length (bp) to be removed (default: %(default)s)',
                        }
                    },
                    {
                        'keys': ['-z', '--gzip'],
                        'properties': {
                            'action': 'store_true',
                            'help': 'gzip the output fastq files',
                        }
                    },
                    WORKDIR_ARG,
                    HELP_ARG,
                    VERSION_ARG,
                ],
        }
}


class EntryPoint:

    root_parser: argparse.ArgumentParser
    variant_filtering_parser: argparse.ArgumentParser
    variant_picking_parser: argparse.ArgumentParser
    vcf2csv_parser: argparse.ArgumentParser
    remove_umi_parser: argparse.ArgumentParser

    def main(self):
        self.set_parsers()
        self.add_arguments()
        self.run()

    def set_parsers(self):
        self.root_parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            formatter_class=argparse.RawTextHelpFormatter,
            add_help=False)

        subparsers = self.root_parser.add_subparsers(
            title='commands',
            dest='mode'
        )

        self.variant_filtering_parser = subparsers.add_parser(
            prog=f'{PROG} {VARIANT_FILTERING}',
            name=VARIANT_FILTERING,
            description=f'{DESCRIPTION} - {VARIANT_FILTERING} mode',
            add_help=False)

        self.variant_picking_parser = subparsers.add_parser(
            prog=f'{PROG} {VARIANT_PICKING}',
            name=VARIANT_PICKING,
            description=f'{DESCRIPTION} - {VARIANT_PICKING} mode',
            add_help=False)

        self.vcf2csv_parser = subparsers.add_parser(
            prog=f'{PROG} {VCF2CSV}',
            name=VCF2CSV,
            description=f'{DESCRIPTION} - {VCF2CSV} mode',
            add_help=False)

        self.remove_umi_parser = subparsers.add_parser(
            prog=f'{PROG} {REMOVE_UMI}',
            name=REMOVE_UMI,
            description=f'{DESCRIPTION} - {REMOVE_UMI} mode',
            add_help=False)

    def add_arguments(self):
        for arg in [HELP_ARG, VERSION_ARG]:
            self.root_parser.add_argument(*arg['keys'], **arg['properties'])

        self.__add(
            parser=self.variant_filtering_parser,
            required_args=MODE_TO_GROUP_TO_ARGS[VARIANT_FILTERING]['Required'],
            optional_args=MODE_TO_GROUP_TO_ARGS[VARIANT_FILTERING]['Optional']
        )

        self.__add(
            parser=self.variant_picking_parser,
            required_args=MODE_TO_GROUP_TO_ARGS[VARIANT_PICKING]['Required'],
            optional_args=MODE_TO_GROUP_TO_ARGS[VARIANT_PICKING]['Optional']
        )

        self.__add(
            parser=self.vcf2csv_parser,
            required_args=MODE_TO_GROUP_TO_ARGS[VCF2CSV]['Required'],
            optional_args=MODE_TO_GROUP_TO_ARGS[VCF2CSV]['Optional']
        )

        self.__add(
            parser=self.remove_umi_parser,
            required_args=MODE_TO_GROUP_TO_ARGS[REMOVE_UMI]['Required'],
            optional_args=MODE_TO_GROUP_TO_ARGS[REMOVE_UMI]['Optional']
        )

    def __add(
            self,
            parser: argparse.ArgumentParser,
            required_args: List[Dict],
            optional_args: List[Dict]):

        group = parser.add_argument_group(f'Required')
        for arg in required_args:
            group.add_argument(*arg['keys'], **arg['properties'])

        group = parser.add_argument_group(f'Optional')
        for arg in optional_args:
            group.add_argument(*arg['keys'], **arg['properties'])

    def run(self):
        args = self.root_parser.parse_args()

        if args.mode is None:
            self.root_parser.print_help()

        elif args.mode == VARIANT_FILTERING:
            print(f'Start running omic {VARIANT_FILTERING} {__VERSION__}\n', flush=True)
            variant_filtering(
                input_vcf=args.input_vcf,
                output_vcf=args.output_vcf,
                variant_flagging_criteria=args.variant_flagging_criteria,
                variant_removal_flags=args.variant_removal_flags,
                only_pass=args.only_pass,
                workdir=args.workdir)

        elif args.mode == VARIANT_PICKING:
            print(f'Start running omic {VARIANT_PICKING} {__VERSION__}\n', flush=True)
            variant_picking(
                ref_fa=args.ref_fa,
                output_vcf=args.output_vcf,
                mutect2=args.mutect2,
                haplotype_caller=args.haplotype_caller,
                muse=args.muse,
                lofreq=args.lofreq,
                varscan=args.varscan,
                vardict=args.vardict,
                somatic_sniper=args.somatic_sniper,
                min_snv_callers=args.min_snv_callers,
                min_indel_callers=args.min_indel_callers,
                workdir=args.workdir)

        elif args.mode == VCF2CSV:
            print(f'Start running omic {VCF2CSV} {__VERSION__}\n', flush=True)
            vcf2csv(
                input_vcf=args.input_vcf,
                output_csv=args.output_csv,
                workdir=args.workdir)

        elif args.mode == REMOVE_UMI:
            print(f'Start running omic {REMOVE_UMI} {__VERSION__}\n', flush=True)
            remove_umi(
                input_fq1=args.input_fq1,
                input_fq2=args.input_fq2,
                output_fq1=args.output_fq1,
                output_fq2=args.output_fq2,
                umi_length=args.umi_length,
                gzip=args.gzip,
                workdir=args.workdir)


if __name__ == '__main__':
    EntryPoint().main()
