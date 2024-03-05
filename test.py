import unittest
import subprocess
from shutil import rmtree


class TestCLI(unittest.TestCase):

    def setUp(self):
        self.workdir = './temp'

    def tearDown(self):
        rmtree(self.workdir)

    def test_variant_filtering(self):
        cmd = f'''python __main__.py variant-filtering \\
--input-vcf ./data/tiny.vcf \\
--output-vcf {self.workdir}/output.vcf \\
--variant-flagging-criteria "LOW_DP: DP<20, HIGH_MQ: MQ>=30" \\
--variant-removal-flags panel_of_normal,LOW_DP \\
--only-pass \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_variant_picking(self):
        cmd = f'''python __main__.py variant-picking \\
--ref-fa ./data/chr9.fa \\
--mutect2 ./data/mutect2.vcf.gz \\
--muse ./data/muse.vcf.gz \\
--lofreq ./data/lofreq.vcf \\
--output-vcf {self.workdir}/output.vcf \\
--min-snv-callers 2 \\
--min-indel-callers 1 \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_vcf2csv(self):
        cmd = f'''python __main__.py vcf2csv \\
--input-vcf ./data/vep.vcf \\
--output-csv {self.workdir}/output.csv \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_remove_umi(self):
        cmd = f'''python __main__.py remove-umi \\
--input-fq1 ./data/tumor.1.fq.gz \\
--input-fq2 ./data/tumor.2.fq.gz \\
--output-fq1 {self.workdir}/output.1.fq.gz \\
--output-fq2 {self.workdir}/output.2.fq.gz \\
--umi-length 7 \\
--gzip \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)
