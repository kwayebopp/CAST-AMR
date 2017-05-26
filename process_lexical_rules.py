#!/usr/bin/python
import sys
from collections import defaultdict
import cPickle
def preprocess(rule_file, result_file):
    rule_count = {}
    with open(rule_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split(' ||| ')
            if ((fields[0], fields[1]) not in rule_count):
                rule_count[(fields[0], fields[1])] = defaultdict(int)

            #diff_part = ' ||| '.join(fields[2:])
            rule_count[(fields[0], fields[1])][fields[2]] += 1
            #rule_count[(fields[0], fields[1])][diff_part] += 1
        f.close()

    with open(result_file, 'w') as wf:
        for (lhs, f_str) in rule_count:
            #str_count_list = sorted(rule_count[(lhs, f_str)].items(), key=lambda elem: (elem[1], len(elem[0].split(' ||| ')[0].split())))
            str_count_list = sorted(rule_count[(lhs, f_str)].items(), key=lambda elem: (elem[1], len(elem[0].split())))

            (selected_e_str, count) = str_count_list[-1]
            wf.write("%s ||| %s ||| %s\n" % (lhs, f_str, selected_e_str))
            if len(str_count_list) >= 2:
                (selected_e_str, count) = str_count_list[-2]
                if count >= 5:
                    wf.write("%s ||| %s ||| %s\n" % (lhs, f_str, selected_e_str))

            #print count
        wf.close()

#Map from a lexicon to its possible candidate fragments
def build_candidates(rule_file, maps_file, max_type):
    rule2cate = defaultdict(set)

    with open(rule_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split(' ||| ')
            #lexicon = fields[1].strip() #Caution here!!
            lexicon = fields[1].strip().lower()

            value = (fields[0].strip(), fields[2].strip())
            if value[0].strip() != 'Nothing':
                try:
                    type = int(value[0][2])
                    if type > max_type:
                        continue
                except:
                    print line
                    continue

            rule2cate[lexicon].add(value)

        f.close()

    wf = open(maps_file, 'wb')
    cPickle.dump(rule2cate, wf)
    wf.close()

def initialize_maps(rule_file, maps_file, max_type):
    rule2cate = {}

    with open(rule_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fields = line.split(' ||| ')
            #lexicon = fields[1].strip() #Caution here!!
            lexicon = fields[1].strip().lower()
            contexts = fields[-1].strip().split()
            for i in xrange(8):
                contexts[i] = contexts[i].lower()
            assert len(contexts) == 12  #4 words, 4 lemmas, 4 POSs
            value = (fields[0].strip(), fields[2].strip())
            if value[0].strip() != 'Nothing':
                try:
                    type = int(value[0][2])
                    if type > max_type:
                        continue
                except:
                    print line
                    continue



            if lexicon not in rule2cate:
                rule2cate[lexicon] = {}

            if (contexts[0], contexts[1], contexts[2], contexts[3]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[0], contexts[1], contexts[2], contexts[3])] = defaultdict(int)
            rule2cate[lexicon][(contexts[0], contexts[1], contexts[2], contexts[3])][value] += 1

            if ((contexts[0], contexts[1]), contexts[2]) not in rule2cate[lexicon]:
                rule2cate[lexicon][((contexts[0], contexts[1]), contexts[2])] = defaultdict(int)
            rule2cate[lexicon][((contexts[0], contexts[1]), contexts[2])][value] += 1

            if (contexts[1], (contexts[2], contexts[3])) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[1], (contexts[2], contexts[3]))] = defaultdict(int)
            rule2cate[lexicon][(contexts[1], (contexts[2], contexts[3]))][value] += 1

            if (contexts[1], contexts[2]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[1], contexts[2])] = defaultdict(int)
            rule2cate[lexicon][(contexts[1], contexts[2])][value] += 1

            if (contexts[0], contexts[1]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[0], contexts[1])] = defaultdict(int)
            rule2cate[lexicon][(contexts[0], contexts[1])][value] += 1

            if (contexts[2], contexts[3]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[2], contexts[3])] = defaultdict(int)
            rule2cate[lexicon][(contexts[2], contexts[3])][value] += 1

            if (contexts[1], '') not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[1], '')] = defaultdict(int)
            rule2cate[lexicon][(contexts[1], '')][value] += 1

            if ('', contexts[2]) not in rule2cate[lexicon]:
                rule2cate[lexicon][('', contexts[2])] = defaultdict(int)
            rule2cate[lexicon][('', contexts[2])][value] += 1

            if '' not in rule2cate[lexicon]:
                rule2cate[lexicon][''] = defaultdict(int)
            rule2cate[lexicon][''][value] += 1

            if (contexts[4], contexts[5], contexts[6], contexts[7]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[4], contexts[5], contexts[6], contexts[7])] = defaultdict(int)
            rule2cate[lexicon][(contexts[4], contexts[5], contexts[6], contexts[7])][value] += 1

            if ((contexts[4], contexts[5]), contexts[6]) not in rule2cate[lexicon]:
                rule2cate[lexicon][((contexts[4], contexts[5]), contexts[6])] = defaultdict(int)
            rule2cate[lexicon][((contexts[4], contexts[5]), contexts[6])][value] += 1

            if (contexts[5], (contexts[6], contexts[7])) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[5], (contexts[6], contexts[7]))] = defaultdict(int)
            rule2cate[lexicon][(contexts[5], (contexts[6], contexts[7]))][value] += 1

            if (contexts[5], contexts[6]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[5], contexts[6])] = defaultdict(int)
            rule2cate[lexicon][(contexts[5], contexts[6])][value] += 1

            if (contexts[4], contexts[5]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[4], contexts[5])] = defaultdict(int)
            rule2cate[lexicon][(contexts[4], contexts[5])][value] += 1

            if (contexts[6], contexts[7]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[6], contexts[7])] = defaultdict(int)
            rule2cate[lexicon][(contexts[6], contexts[7])][value] += 1

            if (contexts[5], '') not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[5], '')] = defaultdict(int)
            rule2cate[lexicon][(contexts[5], '')][value] += 1

            if ('', contexts[6]) not in rule2cate[lexicon]:
                rule2cate[lexicon][('', contexts[6])] = defaultdict(int)
            rule2cate[lexicon][('', contexts[6])][value] += 1

            if (contexts[8], contexts[9], contexts[10], contexts[11]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[8], contexts[9], contexts[10], contexts[11])] = defaultdict(int)
            rule2cate[lexicon][(contexts[8], contexts[9], contexts[10], contexts[11])][value] += 1

            if ((contexts[8], contexts[9]), contexts[10]) not in rule2cate[lexicon]:
                rule2cate[lexicon][((contexts[8], contexts[9]), contexts[10])] = defaultdict(int)
            rule2cate[lexicon][((contexts[8], contexts[9]), contexts[10])][value] += 1

            if (contexts[9], (contexts[10], contexts[11])) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[9], (contexts[10], contexts[11]))] = defaultdict(int)
            rule2cate[lexicon][(contexts[9], (contexts[10], contexts[11]))][value] += 1

            if (contexts[9], contexts[10]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[9], contexts[10])] = defaultdict(int)
            rule2cate[lexicon][(contexts[9], contexts[10])][value] += 1

            if (contexts[8], contexts[9]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[8], contexts[9])] = defaultdict(int)
            rule2cate[lexicon][(contexts[8], contexts[9])][value] += 1

            if (contexts[10], contexts[11]) not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[10], contexts[11])] = defaultdict(int)
            rule2cate[lexicon][(contexts[10], contexts[11])][value] += 1

            if (contexts[9], '') not in rule2cate[lexicon]:
                rule2cate[lexicon][(contexts[9], '')] = defaultdict(int)
            rule2cate[lexicon][(contexts[9], '')][value] += 1

            if ('', contexts[10]) not in rule2cate[lexicon]:
                rule2cate[lexicon][('', contexts[10])] = defaultdict(int)
            rule2cate[lexicon][('', contexts[10])][value] += 1
        f.close()

    wf = open(maps_file, 'wb')
    cPickle.dump(rule2cate, wf)
    wf.close()

if __name__ == '__main__':
    #preprocess(sys.argv[1], sys.argv[2])
    #initialize_maps(sys.argv[1], sys.argv[2], sys.argv[3])
    build_candidates(sys.argv[1], sys.argv[2], sys.argv[3])
