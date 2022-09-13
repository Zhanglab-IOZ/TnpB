import numpy
import pandas
import weblogo

import sys

input_data = pandas.read_table(sys.stdin)
alphabet = weblogo.seq.unambiguous_dna_alphabet
for letter in alphabet.letters():
    if letter not in input_data:
        input_data[letter] = 0
counts = input_data[list(alphabet.letters())].to_numpy()
entropy = input_data["entropy"].to_numpy()
entropy_range = input_data[["lower", "upper"]].to_numpy()

logodata = weblogo.logo.LogoData(
    counts=counts,
    entropy=entropy,
    entropy_interval=entropy_range,
    length=5, alphabet=weblogo.seq.unambiguous_dna_alphabet , weight=numpy.array([1]*5))

logooptions =weblogo. LogoOptions(
    show_fineprint =False,
    color_scheme=weblogo.std_color_schemes["classic"],
    stack_width = 5.4*3,
)
logoformat = weblogo.LogoFormat(logodata, logooptions)
sys.stdout.buffer.write( weblogo.svg_formatter(logodata, logoformat) )
