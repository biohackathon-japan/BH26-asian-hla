"""Order-independent scoring for SpecHLA's unlinked phase blocks.

This preserves the original identity-times-aligned-length objective and API.
It is installed only into an experiment-local copy of SpecHLA.
"""


class Linkage:
    def __init__(self, merged_score_dict):
        self.high_score = 0
        self.support_allele = []
        for allele, scores in sorted(merged_score_dict.items()):
            value = scores[0] * scores[1] + scores[2] * scores[3]
            if value > self.high_score:
                self.high_score = value
                self.support_allele = [[allele] + scores]
            elif value == self.high_score:
                self.support_allele.append([allele] + scores)


def read_blast(self, blast_file):
    """Retain the best HSP per subject, instead of the last BLAST output row."""
    scores = {}
    with open(blast_file) as stream:
        for line in stream:
            row = line.split()
            allele = row[1]
            score = [round(float(row[2]), 2), int(row[3])]
            previous = scores.get(allele)
            key = lambda x: (x[0] * x[1], x[1], x[0])
            if previous is None or key(score) > key(previous):
                scores[allele] = score
    return scores
