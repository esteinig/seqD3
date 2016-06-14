__author__ = 'esteinig'

"""

seqD3
Version 0.1

Multiple Sequence Comparison Tool for interactive visualization in D3
Idea with reference to the BLAST Ring Image Generator (BRIG) by
Alikhan et al. (2012): http://sourceforge.net/projects/brig/

Manual and tutorials for seqD3: https://github.com/esteinig/seqD3.

Eike J. Steinig
Tong Lab, Menzies School of Health Research

eike.steinig@menzies.edu.au

"""

import os
import csv
import json
import numpy
import statistics

from Bio import SeqIO
from Bio import Entrez
from subprocess import call
from html.parser import HTMLParser
from urllib.error import HTTPError
from sklearn.preprocessing import MinMaxScaler

############# Data Module #############

class SeqData:

    """" Data module to download, hold and annotate sequences with Prokka and BioPython. """

    #########################################################################################
    ################                UNDER CONSTRUCTION                #######################
    #########################################################################################

    def __init__(self):

        self.names = []
        self.colors = []

        self.fasta = []
        self.genbank = []

        self.fasta_seqs = []
        self.genbank_seqs = []

    def getGenbank(self, accessions, outdir=os.getcwd(), email="", db='nuccore', rettype="gb", retmode="text"):

        Entrez.email = email

        try:
            seqs = [Entrez.efetch(db=db, rettype=rettype, retmode=retmode, id=id) for id in accessions]
        except HTTPError:
            print("Accession not found.")
        except IOError:
            print("Problem connecting to NCBI.")

    def readFiles(self, files, mode='gb'):

        seqs = [SeqIO.read(open(file, "r"), mode) for file in files]

        if mode == "genbank":
            self.genbank += [file for file in files]
            self.genbank_seqs += seqs
        elif mode == "fasta":
            self.fasta += [file for file in files]
            self.fasta_seqs += seqs

    def reset(self):

        self.fasta_seqs = []
        self.genbank_seqs = []
        self.fasta = []
        self.genbank = []

    def annotateFiles(self, prokka, files=None, path=None):

        prokka.runProkka(files=files, path=path)

        for p in prokka.result_paths:
            gbk_files = [f for f in os.listdir(p) if f.endswith('.gbk')]
            self.readFiles(gbk_files)

        self.readFiles(prokka.files, mode='fasta')

class Prokka:

    def __init__(self):

        self.files = []
        self.path = ''
        self.prefix = ''
        self.locustag = 'PROKKA'
        self.outdir = os.getcwd()

        self.kingdom = 'Bacteria'
        self.genus = ''
        self.species = ''
        self.strain = ''

        self.force = True
        self.addgenes = True
        self.compliant = True
        self.usegenus = True

        self.cpus = 4
        self.evalue = '1e-06'

        self.result_paths = []

    def setOptions(self, outdir=os.getcwd(), force=False, addgenes=False, compliant=False, locustag='PROKKA', cpus=4,
                   evalue='1e-06', kingdom='Bacteria', genus='', species='', strain='', usegenus=False):

        self.locustag = locustag
        self.outdir = outdir

        self.kingdom = kingdom
        self.genus = genus
        self.species = species
        self.strain = strain

        self.force = force
        self.addgenes = addgenes
        self.compliant = compliant
        self.usegenus = usegenus

        self.cpus = str(cpus)
        self.evalue = evalue

    def runProkka(self, files=None, path=None):

        if path is None and files is None:
            raise ValueError('You must provide either a list of files or ' +
                             'the path to a directory containing only the files for annotation with Prokka.')

        if path is None:
            path = ''
            self.files = files
        else:
            self.files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]

        for file in self.files:

            filename = file.split('.')[0]
            print(filename)
            result_path = os.path.join(self.outdir, filename)

            cmd = ['prokka', '--prefix', filename, '--outdir', result_path, '--locustag', self.locustag, '--cpus', self.cpus,
                   '--kingdom', self.kingdom, '--genus', self.genus, '--species', self.species, '--strain', self.strain,
                   '--evalue', self.evalue]

            if self.force:
                cmd.append('--force')
            if self.addgenes:
                cmd.append('--addgenes')
            if self.compliant:
                cmd.append('--compliant')
            if self.usegenus:
                cmd.append('--usegenus')

            cmd.append(os.path.join(path, file))

            call(cmd)

            self.result_paths.append(result_path)




######### Style Handlers #######

# brigD3
# seqD3

############# Generators #############

class RingGenerator:

    """

    Class: Ring Generator

    Initiate with list of rings and set options for visualization. The generator transforms the ring data into
    a list of dictionaries for writing as JSON. Options can be set via setOptions. The main access is the seqD3 method,
    which initiates the visualization class seqD3 and sets its parameters in the string of the script.
    The method then writes HTML file to working directory.

    Attributes:

    self.rings:     list, ring objects
    self.project:   str, name of output files
    self.radius:    int, radius of the ring center
    self.gap:       int, gap between the rings

    self.data:      dict, rings and chords (keys) data for seqD3

    """

    def __init__(self, rings):

        """Initialize generator with list of rings."""

        self.rings = rings
        self.project = 'seqD3'
        self.radius = 300
        self.gap = 10

        self.data = {'rings': [], 'chords': []}

        self.options = {}

        self.ring_animation = """.transition()
                               .delay(function(d, i){ return i * 1; })
                               .duration(5000)
                               .style("opacity", function(d){ return d.opacity; })
                              """

        self.chord_animation = """.transition()
                               .delay(5000)
                               .duration(5000)
                               .style("opacity", function(d){ return d.opacity; })
                               """

        self.text_animation = ".transition().duration(5000).delay(2000)"

    def getData(self):

        """Transform data from rings to data appropriate for D3"""

        lengths = set([ring.length for ring in self.rings])

        if len(lengths) > 1:
            raise ValueError('Ring lengths are not the same.')

        self.options['main_length'] = lengths.pop()

        radius = self.radius

        for ring in self.rings:
            print(radius)
            for seg in ring.data:
                height = seg.pop('height')
                seg['inner'] = float(radius)
                seg['outer'] = float(radius + height)
                self.data['rings'].append(seg)

            radius = radius + ring.height + self.gap

            if isinstance(ring, SequenceRing):
                if ring.chords is not None:
                    self.data['chords'] += ring.chord_data

    def setOptions(self, radius=300, gap=5, project='Test1', title='SCCmec-IV', title_x=40, title_y=40,
                   width=1700, height=1000, animated=True):

        """Set options for circle generator and visualization with D3."""

        self.radius = radius
        self.gap = gap
        self.project = project

        self.options['main_title'] = title
        self.options['main_width'] = width
        self.options['main_height'] = height
        self.options['title_width'] = title_x
        self.options['title_height'] = title_y

        if animated:
            r_ani = self.ring_animation
            c_ani = self.chord_animation
            t_ani = self.text_animation
            s_op = 0
        else:
            r_ani = ''
            c_ani = ''
            t_ani = ''
            s_op = 'function(d, i){ return d.opacity; }'

        self.options['ring_animation'] = r_ani
        self.options['chord_animation'] = c_ani
        self.options['text_animation'] = t_ani
        self.options['start_opacity'] = s_op

    def seqD3(self):

        """Write HTML file for seqD3 to working directory."""

        self.getData()

        print('\nWriting visualization ', self.project + '.html', 'to working directory ...\n')
        viz = seqD3(self.options, self.data)
        viz.setScript()
        viz.writeHTML(self.project)


#  Add manual chord placement.

class ChordGenerator:

    """Generate combinations of comparisons with BLAST and create ChordSets."""

    def __init__(self):

        self.colors = []
        self.mode = 'neighbor'
        self.min_identity = 90
        self.min_length = 2000

        self.chords = []

    def setOptions(self, colors=None, mode='neighbor', min_identity=70, min_length=1000):

        self.mode = mode
        self.colors = colors
        self.min_identity = min_identity
        self.min_length = min_length

    def getChords(self, fastas=None, chord_order=None):

        self.fastas = fastas

        if self.mode == 'neighbor':

            """Neighbor comparisons of .fasta sequences with BLAST"""

            for i in range(len(self.fastas)):

                db = self.fastas[i]
                db_name = db.split('.')[0]

                if i != len(self.fastas)-1:
                    query = self.fastas[i+1]
                    db_idx = i+1
                    src_idx = i+2
                else:
                    query = self.fastas[0]
                    db_idx = i+1
                    src_idx = 1

                query_name = query.split('.')[0]

                blaster = Blaster(db, [query])
                blaster.db = db_name
                blaster.runBLAST()

                chords = ChordSet()
                chords.color = self.colors[i]
                chords.source = src_idx
                chords.target = db_idx
                chords.src_name = query_name
                chords.tar_name = db_name
                chords.min_identity = self.min_identity
                chords.min_length = self.min_length

                chords.readComparison(query_name + 'vs' + db_name)

                self.chords.append(chords)

            if len(fastas) == 2:
                self.chords = [self.chords[0]]

            return self.chords

        elif self.mode == "ordered":

            if chord_order is None:
                # Dictionary Keys: single file to be connected; Dictionary Values: List of files to connect to
                raise ValueError("In ordered mode, please provide a dictionary of chord orders using the given files.")

            for source, targets in chord_order.items():

                db = source
                db_name = db.split(".")[0]

                for target in targets:
                    query = target
                    query_name = query.split(".")[0]
                    db_idx = self.fastas.index(db)+1
                    src_idx = self.fastas.index(query)+1

                    blaster = Blaster(db, [query])
                    blaster.db = db_name
                    blaster.runBLAST()

                    chords = ChordSet()
                    chords.color = self.colors[self.fastas.index(query)]
                    chords.source = src_idx
                    chords.target = db_idx
                    chords.src_name = query_name
                    chords.tar_name = db_name
                    chords.min_identity = self.min_identity
                    chords.min_length = self.min_length

                    chords.readComparison(query_name + 'vs' + db_name)

                    self.chords.append(chords)

            return self.chords

############# Ring Superclass, Ring Subclasses and Readers #############

class Ring:

    """

    Super Class Ring:

    A ring object reads and transforms data into the data shape accesses by the RingGenerator and JS. The standard Ring
    has base attributes colour, name, height and the tooltip style that can be set by the user via setOptions. The
    standard Ring reads a comma-delimited file via readRing containig raw data on each segment (no header) in
    the columns:

    Start Position, End Position, Color, Height and HTML string for Tooltip

    It writes data in the same format via writeRing. Ring objects can also be merged via mergeRings, which adds multiple
    rings' data (list of ring objects) to the current ring object. Attributes of the ring object which called the
    method are retained.

    Subclasses of the ring object have additional attributes pertaining to their function, as well as different readers
    for data files.

    Attributes:

    self.name:      str, name to be shown in tooltips
    self.color:     str, color of ring, hex or name
    self.height:    int, height of ring
    self.tooltip:   obj, tooltip object to set header and text colors

    """

    def __init__(self):

        """Initialize ring object."""

        self.color = 'black'
        self.name = ''
        self.height = 20
        self.opacity = 0.8
        self.tooltip = Tooltip()

        self.positions = []
        self.popups = []
        self.colors = []
        self.heights = []
        self.names = []
        self.opacities = []

        self.data = []

    def mergeRings(self, rings):

        """Merge the current ring with a list of ring objects. Add data of ring objects to current ring."""

        for ring in rings:
            self.data += ring.data

    def setOptions(self, name='Ring', color='black', height=20, opacity=0.8, tooltip=None):

        """Set basic attributes for ring object"""

        self.name = name
        self.color = color
        self.height = height
        self.opacity = opacity

        if tooltip is None:
            self.tooltip = Tooltip()
        else:
            self.tooltip = tooltip

    def getRing(self):

        """Get ring data in dictionary format for Ring Generator and D3."""

        n = len(self.positions)

        print('Generating Ring:', self.name)

        for i in range(n):
            data_dict = {}
            data_dict['start'] = self.positions[i][0]
            data_dict['end'] = self.positions[i][1]
            data_dict['color'] = self.colors[i]
            data_dict['text'] = self.popups[i]
            data_dict['height'] = self.heights[i]
            data_dict['name'] = self.names[i]
            data_dict['opacity'] = self.opacities[i]
            self.data.append(data_dict)

        return self.data

    def writeRing(self, file):

        """Write raw ring data to comma-delimited file."""

        with open(file, 'w') as outfile:
            w = csv.writer(outfile, delimiter=',')

            header = [['Segment Name', 'Start', 'End', 'Color', 'Height', 'Opacity', 'Tooltip Text', 'Tooltip HTML']]
            w.writerows(header)

            d = [[segment['name'], segment['start'], segment['end'], segment['color'], segment['height'],
                  segment['opacity'], strip_tags(segment['text']), segment['text']] for segment in self.data]
            w.writerows(d)

    def readRing(self, file):

        """Read raw ring data from comma-delimited file."""

        self._clear()

        with open(file, 'r') as infile:
            reader = csv.reader(infile)
            header = []
            for row in reader:
                if header:
                    self.names.append(row[0])
                    self.heights.append(float(row[4]))
                    self.colors.append(row[3])
                    self.positions.append([float(row[1]), float(row[2])])
                    self.popups.append(row[7])
                    self.opacities.append(row[5])

                    data = {}
                    data['name'] = row[0]
                    data['start'] = float(row[1])
                    data['end'] = float(row[2])
                    data['color'] = row[3]
                    data['height'] = float(row[4])
                    data['opacity'] = float(row[5])
                    data['text'] = row[7]

                    self.data.append(data)
                else:
                    header = row

    def _clear(self):

        """Clear all ring data."""

        self.heights = []
        self.colors = []
        self.positions = []
        self.popups = []
        self.opacities = []
        self.names = []

        self.data = []

class AnnotationRing(Ring):

    def __init__(self):

        Ring.__init__(self)

        self.feature_types = ['CDS']
        self.extract = {'gene': 'Gene: ', 'product': 'Product: '}

        self.snp_length = 100
        self.intergenic = 'yellow'
        self.synonymous = 'orange'
        self.non_synonymous = 'red'

        self.tooltip = Tooltip()

        self.length = 0

    def readSNP(self, file, n=1):

        """Read SNP data from comma_delimited file (without header): SNP ID, Location, Type, Notes, Ref + N x SNPs"""

        self._clear()

        n += 4  # Index 4 is Reference, n = 0: Reference, n = 1: first Sample

        with open(file, 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                if row[n] != row[4] or n == 4:
                    self.positions.append([int(row[1])-self.snp_length//2, int(row[1])+self.snp_length//2])
                    if row[2] == 'intergenic':
                        self.colors.append(self.intergenic)
                    elif row[2] == 'synonymous':
                        self.colors.append(self.synonymous)
                    elif row[2] == 'non-synonymous':
                        self.colors.append(self.non_synonymous)
                    else:
                        self.colors.append(self.color)

                    self.heights.append(self.height)
                    self.popups.append(self.tooltip.getPopup([('Sequence: ', self.name), ('SNP: ', row[0]),
                                                               ('Location:  ', row[1]), ('Type: ', row[2]),
                                                               ('Note: ', row[3])]))
                    self.opacities.append(self.opacity)
                    self.names.append('')

        self.getRing()

    def readGenbank(self, file, strict=False):

        """Read genbank annotation file and extract relevant features and qualifiers."""

        seq = SeqIO.read(open(file, "r"), "genbank")

        self.length = len(seq)

        features = [feature for feature in seq.features if feature.type in self.feature_types]

        # Fix? #

        if strict:
            clean = []
            for feature in features:
                check = True
                for q in self.extract.keys():
                    if q not in feature.qualifiers:
                        check = False

                if check:
                    clean.append(feature)

            features = clean

        for feature in features:

            self.positions.append([int(feature.location.start), int(feature.location.end)])
            self.colors.append(self.color)
            self.heights.append(self.height)
            self.names.append('')
            self.opacities.append(self.opacity)

            qualifier_texts = []
            for qualifier in self.extract.keys():
                if qualifier in feature.qualifiers:
                    text_tuple = (self.extract[qualifier], ''.join(feature.qualifiers[qualifier]))
                    qualifier_texts.append(text_tuple)

            qualifier_texts.insert(0, ('Location: ', str(feature.location.start) + '-' + str(feature.location.end)))
            qualifier_texts.insert(0, ('Sequence: ', self.name))

            popup = self.tooltip.getPopup(qualifier_texts)

            self.popups.append(popup)

        self.getRing()

    def readSegments(self, segs, length):

        """Manual sequence generation for annotation."""

        self.length = length

        for seq in segs:
            self.positions.append([int(seq['start']), int(seq['end'])])
            self.colors.append(seq['color'])
            self.popups.append(seq['text'])
            self.heights.append(self.height)
            self.names.append(seq['name'])
            self.opacities.append(seq['opacity'])

        self.getRing()

class BlastRing(Ring):

    """Sub-class Blast Ring, for depicting BLAST comparisons against a reference DB"""

    def __init__(self):

        """Initialize super-class ring and attributes for Blast Ring."""

        Ring.__init__(self)

        self.min_identity = 70
        self.min_length = 100
        self.values = []

    def setOptions(self, name='Ring', color='black', height=20, opacity=0.8, tooltip=None, min_identity=0.7,
                   min_length=100):

        self.min_identity = min_identity*100
        self.min_length = min_length

        super(BlastRing, self).setOptions()

    def readComparison(self, file):

        """Reads comparison files from BLAST output (--outfmt 6)"""

        self._clear()

        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                positions = sorted([int(row[8]), int(row[9])])
                if positions[1] - positions[0] >= self.min_length and float(row[2]) >= self.min_identity:
                    self.positions.append(positions)
                    self.values.append(float(row[2]))

        self.colors = [self.color for v in self.values]
        self.heights = [self.height for v in self.values]
        self.opacities = [self.opacity for v in self.values]
        self.names = ['' for v in self.values]

        texts = [[('Sequence: ', self.name), ('BLAST Identity: ', str(v) + '%')] for v in self.values]
        self.popups = [self.tooltip.getPopup(text) for text in texts]

        self.getRing()

class CoverageRing(Ring):

    """Subclass Coverage Ring for depicting coverage matrix across genomes (single or average)."""

    def __init__(self):

        """Initialize super-class ring and attributes for Coverage Ring."""

        Ring.__init__(self)

        self.threshold = 0.96
        self.below = '#E41B17'

    def setOptions(self, name='Ring', color='black', height=20, opacity=0.8, tooltip=None, threshold=0.9,
                   threshold_color='#E41B17'):

        self.threshold = threshold
        self.below = threshold_color

        super(CoverageRing, self).setOptions()

    def readCoverage(self, file, sep='\t', mean=True, rescale=False, n=5):

        """Read coverage matrix from file (with header):
         Segment ID, Start Position, End Position, Value Sample1, Value Sample2..."""

        self._clear()

        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter=sep)

            values = []
            header = []
            texts = []
            for row in reader:
                if header:
                    start = int(row[1])
                    end = int(row[2])
                    if mean:
                        value = statistics.mean([float(v) for v in row[3:]])
                        cov = 'Mean Coverage: '
                    else:
                        value = float(row[n])
                        cov = 'Coverage: '
                    values.append(value)
                    self.positions.append((start, end))
                    color = self.below if value < self.threshold else self.color
                    self.colors.append(color)
                    texts.append([('Sequence: ', self.name), ('Location: ', str(start) + ' - ' + str(end)),
                                 (cov, format(value, ".2f"))])
                    self.opacities.append(self.opacity)
                    self.names.append('')

                else:
                    header = row

        if rescale:
            values = numpy.array(values)
            scaler = MinMaxScaler()
            values = scaler.fit_transform(values)
            values.tolist()

        self.heights = [value*self.height for value in values]
        self.popups = [self.tooltip.getPopup(text) for text in texts]

        self.getRing()

class SequenceRing(Ring):

    def __init__(self, segments, chords=None):

        Ring.__init__(self)

        self.segments = segments
        self.chords = chords
        self.length = 0
        self.radius = 300
        self.gap = 5000

        self.seq_starts = []
        self.chord_data = []

    def setOptions(self, name='Ring', color=None, height=20, opacity=0.8, tooltip=None, chord_radius=300, gap=5000):

        self.radius = chord_radius
        self.gap = gap

        super(SequenceRing, self).setOptions()

    def generateSequence(self):

        self.combineSegments()
        self.addChords()

    def combineSegments(self):

        start = 0
        for seg in self.segments:
            self.positions += [[p+start for p in pos] for pos in seg.positions]
            self.heights += [self.height for h in seg.heights]
            self.popups += seg.popups
            self.colors += seg.colors
            self.names += seg.names
            self.opacities += seg.opacities
            self.seq_starts.append(start)
            start += seg.length + self.gap

        self.length = start

        self.getRing()

    def addChords(self):

        if self.chords is None:
            print('No chords added to Sequence Ring.')
        else:
            for chord in self.chords:
                for i in range(len(chord.source_positions)):

                    src_pos = sorted([int(chord.source_positions[i][0]),
                                     int(chord.source_positions[i][1])])

                    target_pos = sorted([int(chord.target_positions[i][0]),
                                     int(chord.target_positions[i][1])])

                    data = {"source": {"start": src_pos[0]+self.seq_starts[chord.source-1],
                                       "end": src_pos[1]+self.seq_starts[chord.source-1],
                                       "radius": self.radius},
                            "target": {"start": target_pos[0]+self.seq_starts[chord.target-1],
                                       "end": target_pos[1]+self.seq_starts[chord.target-1],
                                       "radius": self.radius},
                            "options": {"color": chord.colors[i],
                                        "text": chord.texts[i]},
                                        "opacity": chord.opacities[i]}

                    self.chord_data.append(data)

############# Chord Class and Reader #############

class ChordSet:

    """Class for holding data for a set of chords."""

    def __init__(self):

        self.source_positions = []
        self.target_positions = []
        self.colors = []
        self.opacities = []
        self.texts = []

        self.tooltip = Tooltip()

        self.opacity = 0.3
        self.color = 'gray'
        self.source = 2  # Query Sequence
        self.target = 1  # DB
        self.src_name = 'Seq 1'
        self.tar_name = 'Seq 2'

        self.min_identity = 70
        self.min_length = 1000

    def readComparison(self, file):

        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                src_pos = sorted([int(row[6]), int(row[7])])
                tar_pos = sorted([int(row[8]), int(row[9])])
                if src_pos[1] - src_pos[0] >= self.min_length or tar_pos[1] - tar_pos[0] >= self.min_length:
                    if float(row[2]) >= self.min_identity:

                        self.source_positions.append(src_pos)
                        self.target_positions.append(tar_pos)
                        self.colors.append(self.color)
                        self.opacities.append(self.opacity)
                        text = [('Comparison: ', self.src_name + ' vs. ' + self.tar_name), ('BLAST Identity: ', str(row[2]) + '%'),
                                ('Location ' + self.src_name + ': ', ' - '.join(str(s) for s in sorted(src_pos))),
                                ('Location ' + self.tar_name + ': ',  ' - '.join(str(s) for s in sorted(tar_pos)))]
                        self.texts.append(self.tooltip.getPopup(text))

############# BLAST #############

class Blaster:

    """

    Class: Blaster

    Convenience module to run BLAST+ (needs to be in $PATH).

    Initialize with a string for the reference genome file (.fasta) and a list of strings for the sequence files
    to be compared (.fasta). Main access is through runBlast. Results attribute holds the string names of the
    output comparison files that can be iterated over to create Blast Rings.

    Attributes:

    self.name_db:       str, name of database to be created from the reference sequence
    self.type:          str, database type, either 'nucl' or 'prot'
    self.mode:          str, type of comparison, either 'blastn' for nucleotide or 'blastp' for protein

    """

    def __init__(self, reference, queries):

        """Initialize blast object with reference and sequence files"""

        self.reference = reference
        self.queries = queries
        self.db = 'ReferenceDB'
        self.type = 'nucl'
        self.mode = 'blastn'

        self.results = []

    def _getDB(self):

        """Run makeblastdb to create reference DB."""

        call(['makeblastdb', '-in', os.path.join(os.getcwd(), self.reference), '-dbtype', self.type, '-out', self.db])
        print('\n')

    def runBLAST(self):

        """"Blast sequence files against reference DB."""

        self._getDB()

        refname = self.reference.split('.')[0]

        for query in self.queries:
            seqname = query.split('.')[0]
            print('Blasting', query, 'against Reference DB ...')
            filename = seqname + 'vs' + refname
            call([self.mode, '-query', os.path.join(os.getcwd(), query), '-db', self.db, '-outfmt', '6', '-out',
                  os.path.join(os.getcwd(), filename)])
            self.results.append(filename)
        print('\n')

############# D3 Visualization #############

class seqD3:

    """Helper class Visualization, holds script for JS D3. Methods to write replace options from Ring Generator in
    script and write the HTML. Initialize with options dict and data from Ring Generator. """

    def __init__(self, options, data):

        self.options = options
        self.data = data

        self.head = """
                    <!DOCTYPE html>

                        <meta charset="UTF-8">
                        <html lang="en">
                        <style>

                        .arcText {
                            fill: #6B6B6B;
                            font-size: 11px;
                            font-family: 'Courgette', sans-serif;
                        }

                        div.tooltip {
                            position: absolute;
                            text-align: left;
                            max-width: 300px;
                            padding: 11px;
                            font-size: 11px;
                            font-family: 'Courgette', sans-serif;
                            background: #FEFCFF;
                            border-radius: 11px;
                            pointer-events: none;
                        }

                        .title {
                            fill: #565051;
                            font-size: 14px;
                            font-family: 'Courgette', sans-serif;
                        }

                        </style>
                        <head>
                            <title></title>
                        </head>
                        <body>
                        <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>
                        """

        self.data1 = """
                        <script type="application/json" id="rings">
                     """

        self.data2 = """
                        </script>
                        <script type="application/json" id="chords">
                     """

        self.script = """
                        </script>
                        <script>
                        var width = main_width
                        var height = main_height
                        var length = main_length
                        var pi = Math.PI
                        var rings = JSON.parse(document.getElementById('rings').innerHTML);
                        var chords = JSON.parse(document.getElementById('chords').innerHTML);

                        // Main Body with Zoom

                        var chart = d3.select("body")
                                        .append("svg")
                                        .attr("id", "chart")
                                        .attr("width", width)
                                        .attr("height", height)
                                        .call(d3.behavior.zoom().on("zoom", function () {
                                            chart.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
                                        }))
                                        .append("g");

                        // Tooltip CSS

                        var div = d3.select("body").append("div")
                                .attr("class", "tooltip")
                                .style("opacity", 0);

                        // Circle Scaler

                        var degreeScale = d3.scale.linear()
                                              .domain([0, length])
                                              .range([0,360]);

                         // Rings

                        var ringShell = chart.append("g")
                                            .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

                        var seqArc = d3.svg.arc(rings)
                                .innerRadius(function(d, i){return d.inner;})
                                .outerRadius(function(d, i){return d.outer;})
                                .startAngle(function(d, i){return degreeScale(d.start) * (pi/180);})
                                .endAngle(function(d, i){return degreeScale(d.end) * (pi/180);})                        ;

                        ringShell.selectAll("path")
                                .data(rings)
                                .enter()
                                .append("path")
                                .style("fill", function(d, i){ return d.color; })
                                .style("opacity", start_opacity)
                                .attr("id", function(d,i) { return "seqArc_"+i; })
                                .attr("d", seqArc)
                                .attr('pointer-events', 'none')
                                .on("mouseover", function(d){
                                        div.transition()
                                            .duration(200)
                                            .style("opacity", .9);
                                        div .html(d.text)
                                            .style("left", (d3.event.pageX + 20) + "px")
                                            .style("top", (d3.event.pageY + 10) + "px");
                                        })
                                .on('mouseout', function(d) {
                                        div.transition()
                                            .duration(200)
                                            .style("opacity", 0)
                                    })

                               ring_animation
                               .transition().attr('pointer-events', 'visible');

                        ringShell.selectAll(".arcText")
                            .data(rings)
                            .enter().append("text")
                            .attr("dy", -13)
                            .style("opacity", 1)
                            .attr("class", "arcText")
                            .append("textPath")
                            .attr("startOffset", '10%')
                            .attr("xlink:href",function(d,i){return "#seqArc_"+i;})
                            .transition().delay(11000)
                            .text(function(d){return d.name; });

                        // Chord Elements

                        var cordShell = chart.append("g")
                                        .attr("class", "chords")
                                        .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

                        var chord = d3.svg.chord(chords)
                                        .radius(function(d, i) {return d.radius;})
                                        .startAngle(function(d, i){return degreeScale(d.start) * (pi/180); })
                                        .endAngle(function(d, i){return degreeScale(d.end) * (pi/180);});

                        cordShell.selectAll("path")
                               .data(chords)
                               .enter()
                               .append("path")
                               .attr("class", "chord")
                               .style("fill", function(d) { return d.options.color; })
                               .style("opacity", start_opacity)
                               .attr("d", chord)
                               .attr('pointer-events', 'none')
                               .on("mouseover", function(){
                                                    d3.select(this).transition().duration(300).style("opacity", 0.8);
                                                })
                               .on("mouseout", function(){
                                            d3.select(this).transition().duration(300).style("opacity",  0.3);
                                        })
                               .on("mousedown", function(d){
                                        div.transition()
                                            .duration(200)
                                            .style("opacity", .9);
                                        div .html(d.options.text)
                                            .style("left", (d3.event.pageX + 20) + "px")
                                            .style("top", (d3.event.pageY + 10) + "px");
                                        })
                                .on('mouseup', function(d) {
                                        div.transition()
                                            .duration(200)
                                            .style("opacity", 0)
                                    })
                               chord_animation
                               .transition().attr('pointer-events', 'visible');

                        // Title

                        titleShell = chart.append("g")
                                         .attr("transform", "translate(" + (width / 100) * title_width + "," + (height / 100) * title_height + ")");

                        titleShell.append("text")
                          .style("opacity", 0)
                          .style("text-anchor", "middle")
                          .style("font-size", "300%")
                          .style("font-weight", "bold")
                          .style("font-family", "times")
                          .attr("class", "title")
                          .text("main_title")
                          text_animation
                          .style("opacity", 1);

                        </script>
                        </body>
                        </html>
                     """

    def setScript(self):

        """Replace placeholder values in script with given options."""

        for placeholder, value in self.options.items():
            self.script = self.script.replace(str(placeholder), str(value))

    def writeHTML(self, project):

        """Write script to HTML."""

        with open(project + '.html', 'w') as outfile:
            outfile.write(self.head)

        with open(project + '.html', 'a') as outfile:
            outfile.write(self.data1)
            json.dump(self.data['rings'], outfile, indent=4, sort_keys=True)
            outfile.write(self.data2)
            json.dump(self.data['chords'], outfile, indent=4, sort_keys=True)
            outfile.write(self.script)


############# Helper Classes #############

class Tooltip:

    """Tooltip class (under construction), will hold more options to customize Tooltips for D3"""

    def __init__(self):

        self.text_color = '#565051'
        self.head_color = '#565051'

    def getPopup(self, text):

        """Converts text - tuple of header and text, i.e. ('Genome:', 'DAR4145') - to HTML string for Tooltip."""

        if len(text) > 0:
            popup = ''
            for i in range(len(text)):
                popup += '<strong>' + '<span style="color:' + self.head_color + '">' + text[i][0] +\
                         '</span>' + '</strong>' + '<span style="color:' + self.text_color + '">' + text[i][1] + ' ' +\
                         '</span>' + '<br>'
        else:
            popup = '-'

        return popup


class MLStripper(HTMLParser):

    def __init__(self):
        super().__init__()

        self.reset()
        self.strict = False
        self.convert_charrefs= True
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        return ''.join(self.fed)

def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    return s.get_data()

