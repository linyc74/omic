import unittest
import subprocess
from os import remove
from shutil import rmtree


class TestCLI(unittest.TestCase):

    def setUp(self):
        self.workdir = './temp'

    def tearDown(self):
        rmtree(self.workdir)
        remove(self.output)

    def test_filtering(self):
        self.output = './output.vcf'
        cmd = f'''python __main__.py filtering \\
--input-vcf ./data/tiny.vcf \\
--output-vcf {self.output} \\
--variant-flagging-criteria "LOW_DP: DP<20, HIGH_MQ: MQ>=30" \\
--variant-removal-flags panel_of_normal,LOW_DP \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_picking(self):
        self.output = './output.vcf'
        cmd = f'''python __main__.py picking \\
--ref-fa ./data/chr9.fa \\
--mutect2 ./data/mutect2.vcf.gz \\
--muse ./data/muse.vcf.gz \\
--lofreq ./data/lofreq.vcf \\
--output-vcf {self.output} \\
--min-snv-callers 2 \\
--min-indel-callers 1 \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_vcf2csv(self):
        self.output = './output.csv'
        cmd = f'''python __main__.py vcf2csv \\
--input-vcf ./data/vep.vcf \\
--output-csv {self.output} \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

