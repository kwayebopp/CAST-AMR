def extract_bigram(tok_list, start, end, win, prefix):
    assert start <= end
    fs = []
    for id in xrange(1, win+1):
        if start - id < 0:
            prev_tok = 'SOS'
        else:
            prev_tok = tok_list[start-id]

        if start - id < -1:
            curr_tok = 'SOS'
        else:
            curr_tok = tok_list[start-id+1]

        curr_f = '%s_lbigram_%d=%s_%s' % (prefix, id, prev_tok, curr_tok)
        fs.append(curr_f)
        #curr_f = '%s_lowerlbigram_%d=%s' % (prefix, id, ('%s_%s' % (prev_tok, curr_tok)).lower())
        #fs.append(curr_f)

    for id in xrange(1, win+1):
        if end + id >= len(tok_list):
            next_tok = 'EOS'
        else:
            next_tok = tok_list[end+id]

        if end + id > len(tok_list):
            curr_tok = 'EOS'
        else:
            curr_tok = tok_list[end+id-1]
        curr_f = '%s_rbigram_%d=%s_%s' % (prefix, id, curr_tok, next_tok)
        fs.append(curr_f)
        #curr_f = '%s_lowerrbigram_%d=%s' % (prefix, id, ('%s_%s' % (curr_tok, next_tok)).lower())
        #fs.append(curr_f)
    return fs

def extract_span(tok_list, start, end, win, prefix, lower=False):
    assert start <= end
    fs = []
    for id in xrange(1, win+1):
        if start - id < 0:
            curr_tok = 'SOS'
        else:
            curr_tok = tok_list[start-id]
        curr_f = '%s_left_%d=%s' % (prefix, id, curr_tok)
        fs.append(curr_f)
        if lower:
            curr_f = '%s_lowerleft_%d=%s' % (prefix, id, curr_tok.lower())
            fs.append(curr_f)

    for id in xrange(1, win+1):
        if end + id >= len(tok_list):
            curr_tok = 'EOS'
        else:
            curr_tok = tok_list[end+id]
        curr_f = '%s_right_%d=%s' % (prefix, id, curr_tok)
        fs.append(curr_f)
        if lower:
            curr_f = '%s_lowerright_%d=%s' % (prefix, id, curr_tok.lower())
            fs.append(curr_f)
    return fs

def extract_curr(tok_list, start, end, prefix, lower=False):
    assert start <= end
    index = 0
    fs = []
    for id in xrange(start, end+1):
        curr_f = '%s_curr_%d=%s' % (prefix, index, tok_list[id])
        fs.append(curr_f)
        if lower:
            curr_f = '%s_lowercurr_%d=%s' % (prefix, index, tok_list[id].lower())
            fs.append(curr_f)
        index += 1
    return fs

def extract_seq_feat(tok_list, start, end, prefix, lower=False):
    fs = []
    fs.append('%s_seq:%s' % (prefix, '_'.join(tok_list[start:end+1])))
    if lower:
        fs.append('%s_lower:%s' % (prefix, ('_'.join(tok_list[start:end+1])).lower()))
    return fs

def suffix(tok):
    fs = []
    for id in xrange(1, 4):
        curr_f = 'suffix_%d=%s' % (id, tok[-id:])
        fs.append(curr_f)
    return fs
