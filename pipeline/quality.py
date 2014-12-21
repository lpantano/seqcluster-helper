from do import run
from tools import splitext_plus

ADAPTER = "TCAGAGTTCTACAGTCCGACGATC"


def run_cuadapt(file_in):
    """remove adapter"""
    adapter = ADAPTER
    output = "%s-clean%s" % (splitext_plus(file_in)[0], ".fastq")
    cmd = ("cutadapt -g {adapter} --discard-untrimmed -O 10 "
            " --minimum-length=25 -o {output} -e 0.2 "
            "{file_in}").format(**locals())
    run(cmd)
    return output
