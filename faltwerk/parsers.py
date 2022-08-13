import os
import numpy as np
import pandas as pd


def get_chunk(stream, separator, read_size=4096):
    """Read from a file chunk by chunk based on a separator substring
    This utility of this function is to avoid reading in the entire contents of a file all at once.
    Instead, you can read in a chunk, process it, then read in the next chunk, and repeat this until
    the EOF.
    Parameters
    ==========
    stream : _io.TextIOWrapper
        A file handle, e.g. stream = open('<path_to_file>', 'r')
    separator : str
        Each value returned will be the string from the last `separator` to the next `separator`
    read_size : int, 4096
        How big should each read size be? Bigger means faster reading, but higher memory usage. This
        has no effect on what is returned, but can greatly influence speed. Default is 4MB.
    References
    ==========
    https://stackoverflow.com/questions/47927039/reading-a-file-until-a-specific-character-in-python
    """

    contents_buffer = ''
    while True:
        chunk = stream.read(read_size)
        if not chunk:
            yield contents_buffer
            break

        contents_buffer += chunk
        while True:
            try:
                part, contents_buffer = contents_buffer.split(separator, 1)
            except ValueError:
                break
            else:
                yield part




class HMMERStandardOutput(object):
    # https://github.com/merenlab/anvio/blob/f9b1de7ee026647a7f7b182567f4f7a8cd61bea3/anvio/parsers/hmmer.py#L25
    """Parse the standard output of HMMER programs (NOTE: currently only works with hmmsearch)
    The main meat of this class is to produce the attributes:
        (1) self.seq_hits
        (2) self.dom_hits
        (3) self.ali_info
    (1) self.seq_hits is a dataframe that looks like this:
        |              query         acc target  query_len        evalue  score  bias  \
        | 0       3Beta_HSD  PF01073.18   1998        282  5.200000e-23   76.2   0.0
        | 1       3Beta_HSD  PF01073.18   1723        282  1.300000e-07   25.7   0.0
        | ...           ...         ...    ...        ...           ...    ...   ...
        | 3128  Voltage_CLC  PF00654.19    320        354  7.200000e-65  214.3  37.1
        | 3129         YkuD  PF03734.13     30        146  1.700000e-14   49.3   0.2
        |       best_dom_evalue  best_dom_score  best_dom_bias  expected_doms  num_doms
        | 0        6.600000e-22            72.6            0.0            2.0         1
        | 1        1.700000e-07            25.3            0.0            1.2         1
        | ...               ...             ...            ...            ...       ...
        | 3128     7.800000e-64           210.9           29.1            2.0         1
        | 3129     3.800000e-14            48.2            0.2            1.7         1
    (2) self.dom_hits is a frame that looks like this:
        |               query         acc target  domain qual  score  bias      c-evalue  \
        | 0       3Beta_HSD  PF01073.18   1998       1    !   72.6   0.0  2.900000e-24
        | 1       3Beta_HSD  PF01073.18   1723       1    !   25.3   0.0  7.300000e-10
        | ...           ...         ...    ...     ...  ...    ...   ...           ...
        | 2896  Voltage_CLC  PF00654.19    320       1    !  210.9  29.1  1.700000e-66
        | 2897         YkuD  PF03734.13     30       1    !   48.2   0.2  8.400000e-17
        |
        |           i-evalue  hmm_start  hmm_stop hmm_bounds  ali_start  ali_stop  \
        | 0     6.600000e-22          1       237         [.          4       243
        | 1     1.700000e-07          1        95         [.          4        92
        | ...            ...        ...       ...        ...        ...       ...
        | 2896  7.800000e-64          3       352         ..         61       390
        | 2897  3.800000e-14          2       146         .]        327       459
        |
        |      ali_bounds  env_start  env_stop env_bounds  mean_post_prob  \
        | 0            ..          4       254         ..            0.74
        | 1            ..          4       148         ..            0.72
        | ...         ...        ...       ...        ...             ...
        | 2896         ..         59       392         ..            0.94
        | 2897         ..        326       459         ..            0.78
        |
        |       match_state_align              comparison_align             sequence_align
        | 0     vvtGggGFlGrrivkeLlrl...  +v+Gg+G++G++ v +L++ ...  LVLGGAGYIGSHAVDQLISK...
        | 1     vvtGggGFlGrrivkeLlrl...  ++ Gg+GFlG++i k L+++...  IIFGGSGFLGQQIAKILVQR...
        | ...                       ...                      ...                      ...
        | 2896  gllagllvkrvapeaagsGi...  g++  +++ r+  + a  G ...  GVVFTYFYTRF-GKNASRGN...
        | 2897  kyivvdlaeqrllvlyengk...  +yi++dl++q++ +++ +gk...  NYIEIDLKDQKM-YCFIDGK...
    If you're confused about the meaning of these columns, please see starting from page 32
    of the HMMER guide http://eddylab.org/software/hmmer/Userguide.pdf. There you will be able
    to with relative ease correlate the column names in these tables to what is described
    meticulously in the tutorial. For example, `best_dom_bias` refers to the the 'bias (best 1
    domain)' column.
    (3) ali_info is a nested dictionary that can be used to access on a per-hit basis which residues
        in a sequence aligned to which residues in the HMM.
    Parameters
    ==========
    hmmer_std_out : str
        Path to output of HMMER.
    context : str, None
        If provided, operations specific to a context will also be carried out. Choose from
        {'interacdome'}
    """

    def __init__(self, hmmer_std_out, context=None):

        self.hmmer_std_out = hmmer_std_out
        self.context = context

        self.set_names()

        self.ali_info = {}

        # This is converted to a dataframe after populating
        self.seq_hits = {
            self.query_col: [],
            self.acc_col: [],
            self.target_col: [],
            self.query_len_col: [],
            'evalue': [],
            'score': [],
            'bias': [],
            'best_dom_evalue': [],
            'best_dom_score': [],
            'best_dom_bias': [],
            'expected_doms': [],
            'num_doms': [],
        }

        self.seq_hits_dtypes = {
            self.query_col: str,
            self.acc_col: str,
            self.target_col: str,
            self.query_len_col: int,
            'evalue': float,
            'score': float,
            'bias': float,
            'best_dom_evalue': float,
            'best_dom_score': float,
            'best_dom_bias': float,
            'expected_doms': float,
            'num_doms': int,
        }

        # This is converted to a dataframe after populating
        self.dom_hits = {
            self.query_col: [],
            self.acc_col: [],
            self.target_col: [],
            'domain': [],
            'qual': [],
            'score': [],
            'bias': [],
            'c-evalue': [],
            'i-evalue': [],
            'hmm_start': [],
            'hmm_stop': [],
            'hmm_bounds': [],
            'ali_start': [],
            'ali_stop': [],
            'ali_bounds': [],
            'env_start': [],
            'env_stop': [],
            'env_bounds': [],
            'mean_post_prob': [],
            'match_state_align': [],
            'comparison_align': [],
            'sequence_align': [],
        }

        self.dom_hits_dtypes = {
            self.query_col: str,
            self.acc_col: str,
            self.target_col: str,
            'domain': int,
            'qual': str,
            'score': float,
            'bias': float,
            'c-evalue': float,
            'i-evalue': float,
            'hmm_start': int,
            'hmm_stop': int,
            'hmm_bounds': str,
            'ali_start': int,
            'ali_stop': int,
            'ali_bounds': str,
            'env_start': int,
            'env_stop': int,
            'env_bounds': str,
            'mean_post_prob': float,
            'match_state_align': str,
            'comparison_align': str,
            'sequence_align': str,
        }

        self.delim_query = '//\n'
        self.delim_seq = '>>'
        self.delim_domain = '=='

        self.load()


    def load(self):

        with open(self.hmmer_std_out) as f:
            for i, query in enumerate(get_chunk(f, separator=self.delim_query, read_size=32768)):

                self.process_query(query)

        self.seq_hits = pd.DataFrame(self.seq_hits).astype(self.seq_hits_dtypes)
        self.dom_hits = pd.DataFrame(self.dom_hits).astype(self.dom_hits_dtypes)

        self.additional_processing()


    def find_line(self, condition):
        for line in self.query_lines[self.line_no:]:
            self.line_no += 1

            if line.startswith('#'):
                continue

            if condition(line):
                return line
        else:
            return False


    def read_lines_until(self, condition, include_last=False, store=True):
        lines = []
        return_value = lines if store else True

        for line in self.query_lines[self.line_no:]:
            self.line_no += 1

            if line.startswith('#'):
                continue

            if condition(line):
                if include_last and store:
                    lines.append(line)

                return lines

            if store:
                lines.append(line)
        else:
            if store:
                return lines
            else:
                return False


    def process_query(self, query):
        if self.delim_seq not in query:
            # This query had no hits
            return

        self.query_lines = query.split('\n')
        self.line_no = 0

        line = self.find_line(lambda line: line.startswith('Query:'))
        line_split = line.split()
        query_name = line_split[1]
        query_len = int(line_split[2][line_split[2].find('=')+1:-1])

        line = self.find_line(lambda line: line.startswith('Accession:'))
        acc = line.split()[1]

        line = self.find_line(lambda line: line.lstrip().startswith('E-value'))
        description_index = line.find('Desc')
        fields = line[:description_index].split() # ignore last 'Description' field

        assert len(fields) == 9, "Please report this on github with your HMMER version"

        self.read_lines_until(lambda line: line.lstrip().startswith('-------'), store=False)
        seq_score_lines = self.read_lines_until(lambda line: line == '')

        num_doms_per_seq = {}

        for seq_score_line in seq_score_lines:
            seq_scores = seq_score_line[:description_index].split()

            self.seq_hits[self.query_col].append(query_name)
            self.seq_hits[self.query_len_col].append(query_len)
            self.seq_hits[self.acc_col].append(acc)
            self.seq_hits['evalue'].append(float(seq_scores[0]))
            self.seq_hits['score'].append(float(seq_scores[1]))
            self.seq_hits['bias'].append(float(seq_scores[2]))
            self.seq_hits['best_dom_evalue'].append(float(seq_scores[3]))
            self.seq_hits['best_dom_score'].append(float(seq_scores[4]))
            self.seq_hits['best_dom_bias'].append(float(seq_scores[5]))
            self.seq_hits['expected_doms'].append(float(seq_scores[6]))
            self.seq_hits['num_doms'].append(int(seq_scores[7]))
            self.seq_hits[self.target_col].append(seq_scores[8])

            num_doms_per_seq[seq_scores[8]] = int(seq_scores[7])

        num_seq_hits = len(seq_score_lines)

        for _ in range(num_seq_hits):
            target_name = self.find_line(lambda line: line.startswith(self.delim_seq)).split()[1]

            if num_doms_per_seq[target_name] == 0:
                continue

            self.line_no += 2
            for __ in range(num_doms_per_seq[target_name]):
                dom_score_summary = self.find_line(lambda line: True).split()

                self.dom_hits[self.query_col].append(query_name)
                self.dom_hits[self.acc_col].append(acc)
                self.dom_hits[self.target_col].append(target_name)
                self.dom_hits['domain'].append(dom_score_summary[0])
                self.dom_hits['qual'].append(dom_score_summary[1])
                self.dom_hits['score'].append(dom_score_summary[2])
                self.dom_hits['bias'].append(dom_score_summary[3])
                self.dom_hits['c-evalue'].append(dom_score_summary[4])
                self.dom_hits['i-evalue'].append(dom_score_summary[5])
                self.dom_hits['hmm_start'].append(dom_score_summary[6])
                self.dom_hits['hmm_stop'].append(dom_score_summary[7])
                self.dom_hits['hmm_bounds'].append(dom_score_summary[8])
                self.dom_hits['ali_start'].append(dom_score_summary[9])
                self.dom_hits['ali_stop'].append(dom_score_summary[10])
                self.dom_hits['ali_bounds'].append(dom_score_summary[11])
                self.dom_hits['env_start'].append(dom_score_summary[12])
                self.dom_hits['env_stop'].append(dom_score_summary[13])
                self.dom_hits['env_bounds'].append(dom_score_summary[14])
                self.dom_hits['mean_post_prob'].append(dom_score_summary[15])

            for __ in range(num_doms_per_seq[target_name]):
                self.find_line(lambda line: line.lstrip().startswith(self.delim_domain))

                if __ == num_doms_per_seq[target_name] - 1:
                    if _ == num_seq_hits - 1:
                        # This is the last alignment in the summary_info. Go to end of string
                        ali_lines = self.read_lines_until(lambda line: False)
                    else:
                        # This is the last alignment in the sequence. Go to next sequence delimiter
                        ali_lines = self.read_lines_until(lambda line: line.lstrip().startswith(self.delim_seq))
                        self.line_no -= 1
                else:
                    ali_lines = self.read_lines_until(lambda line: line.lstrip().startswith(self.delim_domain))
                    self.line_no -= 1

                consensus = []
                match = []
                target = []
                line_index = 0
                while True:
                    if line_index >= len(ali_lines):
                        break

                    line = ali_lines[line_index]

                    if not line.lstrip().startswith(query_name + ' '):
                        line_index += 1
                        continue

                    cons_seq_fragment = line.split()[2]
                    frag_len = len(cons_seq_fragment)
                    ali_index = line.find(cons_seq_fragment)

                    consensus.append(cons_seq_fragment)
                    match.append(ali_lines[line_index + 1][ali_index: ali_index + frag_len])
                    target.append(ali_lines[line_index + 2][ali_index: ali_index + frag_len])

                    line_index += 2

                self.dom_hits['match_state_align'].append(''.join(consensus))
                self.dom_hits['comparison_align'].append(''.join(match))
                self.dom_hits['sequence_align'].append(''.join(target))


    def set_names(self):
        """Set the column names depending on self.context"""

        if self.context is None:
            self.query_col = 'query'
            self.acc_col = 'acc'
            self.query_len_col = 'query_len'
            self.target_col = 'target'

        elif self.context == 'interacdome':
            self.query_col = 'pfam_name'
            self.acc_col = 'pfam_id'
            self.query_len_col = 'pfam_len'
            self.target_col = 'corresponding_gene_call'


    def additional_processing(self):
        """Further process raw data"""

        if self.context is None:
            self.get_ali_info()

        elif self.context == 'interacdome':
            self.seq_hits['corresponding_gene_call'] = self.seq_hits['corresponding_gene_call'].astype(int)
            self.dom_hits['corresponding_gene_call'] = self.dom_hits['corresponding_gene_call'].astype(int)

            if self.dom_hits.empty:
                self.dom_hits['version'] = []
            else:
                self.dom_hits[['pfam_id', 'version']] = self.dom_hits['pfam_id'].str.split('.', n=1, expand=True)

            if self.seq_hits.empty:
                self.seq_hits['version'] = []
            else:
                self.seq_hits[['pfam_id', 'version']] = self.seq_hits['pfam_id'].str.split('.', n=1, expand=True)

            # For convenience this is done after pfam_id has been split
            self.get_ali_info()


    def get_ali_info(self):
        """Creates self.ali_info. See class docstring for description
        Notes
        =====
        - This function is very slow.
        - EDIT: This function is not _that_ slow
        """

        if self.dom_hits.empty:
            return

        unique_targets = self.dom_hits[self.target_col].nunique()

        gap_chars = {'-', '.'}

        processed = 0
        for target, subset in self.dom_hits.groupby(self.target_col):

            self.ali_info[target] = {}

            for acc, subsubset in subset.groupby(self.acc_col):
                for i, row in subsubset.iterrows():
                    seq_positions, seq_chars, hmm_positions, hmm_chars, comparison_chars = [], [], [], [], []

                    seq_pos, hmm_pos = row['ali_start'], row['hmm_start']
                    sequence, match_state, comparison = row['sequence_align'], row['match_state_align'], row['comparison_align']

                    assert len(sequence) == len(match_state)

                    for i in range(len(sequence)):
                        seq_char, hmm_char, comparison_char = sequence[i], match_state[i], comparison[i]
                        if (seq_char not in gap_chars) and (hmm_char not in gap_chars):
                            # there is alignment (non-gap characters)
                            seq_positions.append(seq_pos)
                            seq_chars.append(seq_char)
                            hmm_positions.append(hmm_pos)
                            hmm_chars.append(hmm_char.upper())
                            comparison_chars.append(comparison_char.upper())
                            seq_pos += 1
                            hmm_pos += 1
                        elif (seq_char in gap_chars) and (hmm_char not in gap_chars):
                            # gap in seq
                            hmm_pos += 1
                        elif (seq_char not in gap_chars) and (hmm_char in gap_chars):
                            # gap in match state
                            seq_pos += 1
                        else:
                            # this happens with 0 probability
                            pass

                    # The HMM state and sequence positions are 1-indexed. We subtract by 1 to make them zero-indexed
                    self.ali_info[target][(acc, row['domain'])] = pd.DataFrame({
                        'seq': seq_chars,
                        'hmm': hmm_chars,
                        'comparison': comparison_chars,
                        'seq_positions': np.array(seq_positions) - 1,
                        'hmm_positions': np.array(hmm_positions) - 1,
                    })

            processed += 1
