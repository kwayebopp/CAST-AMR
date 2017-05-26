import dill as pickle
import nltk
from nltk.model import *
from nltk.corpus import inaugural, brown

print 'building vocabulary...'
vocab = build_vocabulary(2, brown.words()[0:], inaugural.words()[0:])
print 'vocabulary complete!'


def ct_ngrams(order):
    return count_ngrams(order, vocab, brown.sents()[0:], inaugural.sents()[0:])

def ngram_model(order):
    return LaplaceNgramModel(ct_ngrams(order))
print 'building lang models...'

bg_model = ngram_model(2)
with open('bigram.lm', 'w+b') as f:
    print 'pickling...'
    pickle.dump(bg_model, f, protocol=pickle.HIGHEST_PROTOCOL)
print '33% complete!'
tg_model = ngram_model(3)
with open('trigram.lm', 'w+b') as f:
    print 'pickling...'
    pickle.dump(tg_model, f, protocol=pickle.HIGHEST_PROTOCOL)
print '66% complete!'
qg_model = ngram_model(4)
with open('4gram.lm', 'w+b') as f:
    print 'pickling...'
    pickle.dump(qg_model, f, protocol=pickle.HIGHEST_PROTOCOL)
print 'lang models complete!'