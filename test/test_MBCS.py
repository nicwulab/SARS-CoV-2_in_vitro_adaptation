import os
import filecmp

TEST_PATH = os.path.dirname(os.path.abspath(__file__))
CODE_PATH = TEST_PATH + '/../codes'

def test_extract_MBCS():
    outfile = '{TEST_PATH}/data/test.tsv'.format(TEST_PATH = TEST_PATH)
    cmd = 'python \'{CODE_PATH}/extract_MBCS.py\' \'{TEST_PATH}/data/test.bam\' \'{OUT}\' test.png'.format(CODE_PATH=CODE_PATH, TEST_PATH=TEST_PATH, OUT=outfile)
    print(cmd)
    os.system(cmd)
    assert(filecmp.cmp(outfile, TEST_PATH + '/data/mbcs.out'))
