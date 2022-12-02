import os
import shutil
import unittest
import subprocess


class TestCLI(unittest.TestCase):

    def setUp(self):
        self.workdir = './temp'
        self.output_vcf = './output.vcf'

    def tearDown(self):
        shutil.rmtree(self.workdir)
        os.remove(self.output_vcf)

    def test_filtering(self):
        cmd = f'''python __main__.py filtering \\
--input-vcf ./data/tiny.vcf \\
--output-vcf {self.output_vcf} \\
--variant-flagging-criteria "LOW_DP: DP<20, HIGH_MQ: MQ>=30" \\
--variant-removal-flags panel_of_normal,LOW_DP \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)

    def test_picking(self):
        cmd = f'''python __main__.py picking \\
--ref-fa ./data/chr9.fa \\
--mutect2 ./data/mutect2.vcf \\
--muse ./data/muse.vcf \\
--lofreq ./data/lofreq.vcf \\
--output-vcf {self.output_vcf} \\
--min-snv-callers 2 \\
--min-indel-callers 1 \\
--workdir {self.workdir}'''
        subprocess.check_call(cmd, shell=True)



