from os import popen
from glob import glob
from subprocess import call
from argparse import ArgumentParser

parser = ArgumentParser(
    description='Script to automate running over all input files')
parser.add_argument('-l', '--lepton', action='store',
                    default='e', dest='lepton',
                    help='Which channel? [e, m]'
                    )
parser.add_argument('-e', '--exe', action='store',
                    default='ZTT_XSection_mu.exe', dest='exe',
                    help='name of executable to run'
                    )
parser.add_argument('-i', '--input-dir', action='store',
                    default='/store/user/tmitchel/DAS2019-ZTT-Long/', dest='input_dir',
                    help='path to input files'
                    )
parser.add_argument('-o', '--output-suffix', action='store',
                    default='test_v1', dest='output_suffix',
                    help='suffix to add to output file name'
                    )
args = parser.parse_args()

fileList = [ifile for ifile in filter(None, popen(
    'xrdfs root://cmseos.fnal.gov/ ls {}'.format(args.input_dir)).read().split('\n')) if '.root' in ifile]

for ifile in fileList:
    fname = ifile.split('/')[-1].split('.root')[0]
    call('./{} {}_{}.root root://cmseos.fnal.gov/{}'.format(args.exe,
                                                            fname, args.output_suffix, args.input_dir), shell=True)
