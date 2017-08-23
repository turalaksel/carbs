import argparse
from os import path
import sys
import cadnano
import numpy as np
from cadnano.document import Document

app = cadnano.app()
doc = app.document = Document()
doc.readFile('6hb.json');

part = doc.activePart()

vh0 = np.concatenate(part.getCoordinates(0))



oligos = part.oligos()

oligos_sorted_by_length = sorted(oligos, key=lambda x: x.length(), reverse=True)
longest_oligo = oligos_sorted_by_length[0]
staple_oligos = oligos_sorted_by_length[1:]


for staple in staple_oligos:
    # print(staple)
    for strand in staple.strand5p().generator3pStrand():
        print(strand)


staple
