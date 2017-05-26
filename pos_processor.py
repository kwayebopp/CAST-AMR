#!/usr/bin/python
import nltk
import argparse
def readPOSs(pos_file):
    pos_seqs = []
    with open(pos_file, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue
            word_tags = [('/'.join(word_tag.split('/')[:-1]), word_tag.split('/')[-1]) for word_tag in line.strip().split()]
            pos_seqs.append(word_tags)
        f.close()
    return pos_seqs

def pos_tagging(tok_file, pos_file, separator='##**##'):
    with open(tok_file, 'r') as f:
        with open(pos_file, 'w') as wf:
            for line in f:
                toks = line.decode('utf-8').strip().split()
                seq = []
                for (tok, pos) in nltk.pos_tag(toks):
                    assert ' ' not in tok and ' ' not in pos
                    seq.append(pos)
                assert len(seq) == len(toks)
                print >>wf, ' '.join(seq)
            wf.close()
            f.close()

if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--tok", type=str, help="the original tok file", required=True)
    argparser.add_argument("--output", type=str, help="the output file", required=True)

    args = argparser.parse_args()
    pos_tagging(args.tok, args.output)
