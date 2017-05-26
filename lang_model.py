import nltk
import sys

from nltk.model import *

sents = [line for line in  open('./dev/amr.lin', 'r')]
words = []
for sent in sents: 
    words.append("<s>")
    for word in sent.split():
        if word != 'multi-sentence':
            words.append(word)
    words.append("</s>")


def bigram_model():
    return MLENgramModel(count_ngrams(4, build_vocabulary(2, words[0:]), sents[0:]))


