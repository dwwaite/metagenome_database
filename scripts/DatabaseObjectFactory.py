'''

'''
from collections import Counter

class Contig():

    def __init__(self, contig_name, contig_sequence):

        self._contig_name = contig_name
        self._contig_sequence = contig_sequence
        self._contig_length = len(contig_sequence)
        self._contig_gc = self._calc_gc(contig_sequence)

    def _calc_gc(self, sequence):

        nuclotide_counter = Counter(sequence)

        gc_frac = float( nuclotide_counter['G'] + nuclotide_counter['C'] ) / sum( nuclotide_counter.values() )
        return gc_frac * 100
    @property
    def name(self):
        return self._contig_name

    @property
    def sequence(self):
        return self._contig_sequence

    @property
    def length(self):
        return self._contig_length

    @property
    def gc(self):
        return self._contig_gc

class Gene():

    def __init__(self, gene_name, contig_name, prediction_tool, sequence_nt, sequence_aa, start, stop, orientation):
        self._gene_name = gene_name
        self._contig_name = contig_name
        self._prediction_tool = prediction_tool
        self._sequence_nt = sequence_nt
        self._sequence_aa = sequence_aa
        self._start = start
        self._stop = stop
        self._orientation = orientation

    @property
    def name(self):
        return self._gene_name

    @property
    def contig_name(self):
        return self._contig_name

    @property
    def prediction_tool(self):
        return self._prediction_tool

    @property
    def sequence_nt(self):
        return self._sequence_nt

    @property
    def length_nt(self):
        return len(self._sequence_nt)

    @property
    def sequence_aa(self):
        return self._sequence_aa

    @property
    def length_aa(self):
        return len(self._sequence_aa)

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def orientation(self):
        return self._orientation

class Annotation():

    def __init__(self, gene_name, method, annotation_database, accession, identity, coverage, evalue, description):
        self._gene_name = gene_name
        self._method = method
        self._annotation_database = annotation_database
        self._accession = accession
        self._identity = identity
        self._coverage = coverage
        self._evalue = evalue
        self._description = description
        self._is_current = True

    @property
    def gene_name(self):
        return self._gene_name

    @property
    def method(self):
        return self._method

    @property
    def annotation_database(self):
        return self._annotation_database

    @property
    def accession(self):
        return self._accession

    @property
    def identity(self):
        return self._identity

    @property
    def coverage(self):
        return self._coverage

    @property
    def evalue(self):
        return self._evalue

    @property
    def description(self):
        return self._description

    @property
    def is_current(self):
        return self._is_current

class Mapping():

    def __init__(self, target_name, sample_name, reads_mapped):
        self._target_name = target_name
        self._sample_name = sample_name
        self._reads_mapped = reads_mapped

    @property
    def target_name(self):
        return self._target_name

    @property
    def sample_name(self):
        return self._sample_name

    @property
    def reads_mapped(self):
        return self._reads_mapped

class Bin():

    def __init__(self, bin_name, taxonomy, completeness, contamination, strain_heterogeneity, version_name, date_created):

        self._bin_name = bin_name
        self._taxonomy = taxonomy
        self._completeness = completeness
        self._contamination = contamination
        self._strain_heterogeneity = strain_heterogeneity
        self._version_name = version_name
        self._date_created = date_created

    @property
    def name(self):
        return self._bin_name

    @property
    def taxonomy(self):
        return self._taxonomy

    @property
    def completeness(self):
        return self._completeness

    @property
    def contamination(self):
        return self._contamination

    @property
    def strain_heterogeneity(self):
        return self._strain_heterogeneity

    @property
    def version_name(self):
        return self._version_name

    @property
    def date_created(self):
        return self._date_created