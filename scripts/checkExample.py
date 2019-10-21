#!/usr/bin/env python3

"""
A module which contains all the non-regression tests and unitary tests for intronSeeker's 
python scripts which cannot be wrapped in scripts' docstrings and performed correctly by doctest module.
To run the all the tests use pytest in scripts directory or run this scripts :

:Example:

  ~/intronSeeker/scripts $ py.test --doctest-modules *
  
  OR 
  
  ~$  intronSeeker --checkExample

.. seealso:: checkInstall.py, seekerToolbox.py, dataSimulation.py, readsMapping.py, intronSearch.py
"""

# Metadata
__author__ = 'Lasguignes Emilien - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2019 INRA'
__license__ = 'GNU General Public License'
__version__ = '1.0'
__email__ = 'support.bioinfo.genotoul@inra.fr'
__status__ = 'dev'

from dataSimulation import * # TODO : A remplacer par dataSimulation quand intégré
import os
import tempfile as tmp
import configparser
import subprocess as sp
import pandas as pd
import re
from pprint import pprint

def test_read_gtf() :
    gtf, transcripts = read_gtf(os.path.abspath(os.path.dirname(__file__))+"/../data/test_GBS.gtf")
    grep = sp.Popen(["grep","-v","^#",os.path.abspath(os.path.dirname(__file__))+"/../data/test_GBS.gtf"], stdout=sp.PIPE) ;
    wcline = sp.Popen(["wc", "-l"], stdin=grep.stdout, stdout=sp.PIPE) ;
    grep.stdout.close() ; 
    sdo = wcline.communicate()[0] ;
    assert len(gtf) == int(re.search(r'(\d+)',str(sdo)).group(1))
    grep2 = sp.Popen(["grep", 'transcrip.* "protein_coding";',os.path.abspath(os.path.dirname(__file__))+"/../data/test_GBS.gtf"], stdout=sp.PIPE) ;
    awk = sp.Popen(["awk","-F",'";',"{print $3}"],stdin=grep2.stdout, stdout=sp.PIPE) ;
    sort = sp.Popen(["sort","-u"],stdin=awk.stdout, stdout=sp.PIPE) ;
    wcline2 = sp.Popen(["wc","-l"],stdin=sort.stdout, stdout=sp.PIPE) ;
    sdo = wcline2.communicate()[0]
    assert len(transcripts) == int(re.search(r'(\d+)',str(sdo)).group(1))
    
    

def test_make_density_law() :
    law = make_density_law()
    
    config = configparser.RawConfigParser()
    config.read(os.path.abspath(os.path.dirname(__file__))+"/../config/intronSeeker.properties")
    assert set(law) == {int(k) for k in config["Density"].keys()} # Check if the classes are correct
    for classe in config["Density"] :
        assert law.count(int(classe)) == int(config["Density"][classe]) # Check if the effective for each class is correct

def test_choose_transcripts() :
    transcripts = ['3R5.1a.1', '3R5.1b.1', 'B0348.1.1', 'B0348.10.1', 'B0348.2a.1', 'B0348.2b.1', 'B0348.4a.1']
    law = make_density_law() 
    choosen = choose_transcripts(transcripts,3)
    assert len(choosen) == 3 
    assert set(choosen).issubset(set(transcripts))
    assert set(choosen.values()).issubset(set(law))
    choosen2 = choose_transcripts(transcripts,100)
    assert set(choosen2) == set(transcripts)
    assert set(choosen2.values()).issubset(set(law))

def test_parse_gtf_content() :
    gtf,transcripts = read_gtf(os.path.abspath(os.path.dirname(__file__))+"/../data/test_GBS.gtf")
    choosen = choose_transcripts(transcripts,1000)
    reference, library, control = parse_gtf_content(gtf,dict(choosen),False)
    assert { f[-1]["transcript_id"] for f in reference } == set(choosen)
    assert { f[-1]["transcript_id"] for f in library } == set(choosen)
    for choosen_one in choosen :
        ref_one = [f for f in reference if f[-1]["transcript_id"] == choosen_one and f[2]=="exon"]  
        lib_one = [f for f in library if f[-1]["transcript_id"] == choosen_one and f[2]=="exon"]
        gtf_one = [f for f in gtf if ("transcript_id" in f[-1] and f[-1]["transcript_id"] == choosen_one and f[2]=="exon") ] 
        # Check if all the exons are correctly gathered
        # by checking the correct number of exons and if there is dupicated exons
        # Check also, according to the real class of the transcript,
        # if the modified transcript (with spliced exon or retained intron) is correctly 
        # added to result lists (respectively library or reference).
        if int(ref_one[0][-1]["classe"]) == 0 :
            assert len(ref_one) == len(gtf_one)
            assert { e[-1]["exon_id"] for e in ref_one} == { e[-1]["exon_id"] for e in gtf_one}
            assert len(lib_one) == len(gtf_one)
            assert { e[-1]["exon_id"] for e in lib_one} == { e[-1]["exon_id"] for e in gtf_one}
        elif int(ref_one[0][-1]["classe"]) == -1  :
            assert len(ref_one) == len(gtf_one)
            assert { e[-1]["exon_id"] for e in ref_one} == { e[-1]["exon_id"] for e in gtf_one}
            assert len(lib_one) == len(gtf_one)-1
            assert { e[-1]["exon_id"] for e in lib_one}.issubset({ e[-1]["exon_id"] for e in gtf_one})
        elif int(ref_one[0][-1]["classe"]) > 0 :
            assert len(ref_one) == len(gtf_one)+int(ref_one[0][-1]["classe"])
            assert { e[-1]["exon_id"] if "exon_id" in e[-1] else e[-1]["intron_id"] for e in ref_one}.issuperset({ e[-1]["exon_id"] for e in gtf_one})
            assert len(lib_one) == len(gtf_one)
            assert { e[-1]["exon_id"] for e in lib_one} == { e[-1]["exon_id"] for e in gtf_one}
    # Check if control list generation is correct
    real_class = {f[-1]["transcript_id"] : int(f[-1]["classe"]) for f in reference}
    nb_features = sum([abs(c) for c in real_class.values()])
    assert len(control) == nb_features
    assert len({(f[-1]["transcript_id"],f[-1]["exon_id"]) if "exon_id" in f[-1] else (f[-1]["transcript_id"],f[-1]["intron_id"]) for f in control}) == nb_features

def test_transcript_df() :
    # Check if it generates a correct DataFrame for a transcript on strand + 
    transcript = [ 
        ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1'}], 
        ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2'}]
        ]
    df = pd.DataFrame([
            ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1'}, True],
            ['III', 'WormBase', 'intron', '13782935', '13783360', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'intron_id': '3R5.2.e1+3R5.2.e2'}, False],
            ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2'}, True]
            ],
        columns=["ref","DB","feature","start","end","score","strand","frame","misc_attr","in_transcript"]
        )
    pd.util.testing.assert_frame_equal(transcript_df(transcript),df)
    
    # Check if it generates a correct DataFrame for a transcript on strand -
    transcript = [
    ['V', 'WormBase', 'exon', '13221', '13369', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '1', 'exon_id': 'B0348.9.e1'}],
    ['V', 'WormBase', 'exon', '12525', '12595', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '2', 'exon_id': 'B0348.9.e2'}]
    ]
    df = pd.DataFrame([
            ['V', 'WormBase', 'exon', '13221', '13369', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '1', 'exon_id': 'B0348.9.e1'}, True],
            ['V', 'WormBase', 'intron', '12596', '13220', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'intron_id': 'B0348.9.e1+B0348.9.e2'}, False],
            ['V', 'WormBase', 'exon', '12525', '12595', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '2', 'exon_id': 'B0348.9.e2'}, True]
            ],
        columns=["ref","DB","feature","start","end","score","strand","frame","misc_attr","in_transcript"]
        )
    pd.util.testing.assert_frame_equal(transcript_df(transcript),df)
    
def test_construct_new_transcript() :
    transcript_pos = [ 
        ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1'}], 
        ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2'}]
        ]
    transcript_neg = [
        ['V', 'WormBase', 'exon', '13452', '13809', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '1', 'exon_id': 'B0348.9.e1'}],
        ['V', 'WormBase', 'exon', '13221', '13369', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '2', 'exon_id': 'B0348.9.e2'}],
        ['V', 'WormBase', 'exon', '12525', '12595', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '3', 'exon_id': 'B0348.9.e3'}]
        ]
    
    # Check transcript class 0 
    reference,library,feature = construct_new_transcript(transcript_pos,0)
    assert reference == [ 
            ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1','classe': '0'}], 
            ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2','classe': '0'}]
            ]
    assert library == [ 
            ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1','classe': '0'}], 
            ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2','classe': '0'}]
            ]
    assert not feature
    
    # Check transcript class -1 on strand +
    reference,library,feature = construct_new_transcript(transcript_pos,-1)
    assert reference == [ 
            ['III', 'WormBase', 'exon', '13782587', '13782934', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1','classe': '-1'}], 
            ['III', 'WormBase', 'exon', '13783361', '13783459', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2','classe': '-1'}]
            ]
    assert len(library) == len(transcript_pos)-1
    assert library[0][-1]['classe'] == '-1'
    del library[0][-1]['classe']
    assert library[0] in transcript_pos
    assert len(feature) == 1
    potential_features = {
        '3R5.2.e1' : [['3R5.2', 'WormBase', 'spliced_exon', '1', '348', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '1', 'exon_id': '3R5.2.e1', 'classe': '-1'}]], 
        '3R5.2.e2': [['3R5.2', 'WormBase', 'spliced_exon', '349', '447', '.', '+', '.', {'gene_id': 'WBGene00007066', 'transcript_id': '3R5.2', 'exon_number': '2', 'exon_id': '3R5.2.e2', 'classe': '-1'}]]
        }
    assert feature == potential_features[feature[0][-1]['exon_id']]
    
    # Check transcript class -1 on strand -
    reference,library,feature = construct_new_transcript(transcript_neg,-1)
    assert reference == [
            ['V', 'WormBase', 'exon', '13452', '13809', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '1', 'exon_id': 'B0348.9.e1','classe': '-1'}],
            ['V', 'WormBase', 'exon', '13221', '13369', '.', '-', '.', {'gene_id': 'WBGene00235314','transcript_id': 'B0348.9', 'exon_number': '2', 'exon_id': 'B0348.9.e2','classe': '-1'}],
            ['V', 'WormBase', 'exon', '12525', '12595', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '3', 'exon_id': 'B0348.9.e3','classe': '-1'}]
            ]
    assert len(library) == len(transcript_neg)-1
    assert library[0][-1]['classe'] == '-1'
    del library[0][-1]['classe']
    assert library[0] in transcript_neg
    assert len(feature) == 1
    potential_features = {
        'B0348.9.e1' : [['B0348.9', 'WormBase', 'spliced_exon', '1', '358', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '1', 'exon_id': 'B0348.9.e1', 'classe': '-1'}]],
        'B0348.9.e2' : [['B0348.9', 'WormBase', 'spliced_exon', '359', '507', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '2', 'exon_id': 'B0348.9.e2', 'classe': '-1'}]],
        'B0348.9.e3': [['B0348.9', 'WormBase', 'spliced_exon', '508', '578', '.', '-', '.', {'gene_id': 'WBGene00235314', 'transcript_id': 'B0348.9', 'exon_number': '3', 'exon_id': 'B0348.9.e3', 'classe': '-1'}]]
        }
    assert feature == potential_features[feature[0][-1]['exon_id']]

# ~ import os ;
# ~ import sys ;
# ~ import shutil ;

# ~ def _non_regression_tests() :
	
	# ~ data_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__))+"/../data") ;
	# ~ print(data_path)
	
	# ~ try :
		# ~ os.mkdir("test_dir") ;
	# ~ except FileExistsError :
		# ~ pass
	# ~ finally :
		# ~ os.chdir("test_dir") ;
	
	# ~ split_command = "intronStalker split -i {alignment} -r {reference} -o {output}".format(
		# ~ alignment = data_path + "/Alignment_Example.Aligned.sortedByCoord.out.bam",
		# ~ reference = data_path + "/Reference_Example.fa",
		# ~ output = "checkExample_split" 
		# ~ );

	# ~ print("### Splicing events candidates searching test ###") ;
	# ~ print("Execution of the command :") ;
	# ~ print( split_command +"\n") ;
	# ~ os.system( split_command ) ;
	
	
	#~ print("### Genome-based data simulation test ###") ;
	#~ print("Execution of the command :") ;
	#~ print("intronStalker annoToData -i {genome_gff} -f {genome_fasta} -o regression_test_GBS".format(
			#~ genome_gff=os.path.abspath(os.path.dirname(sys.argv[0]) + "/../data/Cele/GCF_000002985.6_WBcel235_genomic.gff"),
			#~ genome_fasta=os.path.abspath(os.path.dirname(sys.argv[0]) + "/../data/Cele/Cele_whole_genome.fna")
			#~ )) ;
	#~ os.system("intronStalker annoToData -i {genome_gff} -f {genome_fasta} -o regression_test_GBS".format(
			#~ genome_gff=os.path.abspath(os.path.dirname(sys.argv[0]) + "/../data/Cele/GCF_000002985.6_WBcel235_genomic.gff"),
			#~ genome_fasta=os.path.abspath(os.path.dirname(sys.argv[0]) + "/../data/Cele/Cele_whole_genome.fna")
			#~ )) ;
	
	# ~ os.chdir("..") ;
	# ~ shutil.rmtree("test_dir") ;

#~ def _check_command_output(dir_name : str, nb_file : int) :
	#~ try :
		#~ if dir_name not in os.listdir(".") :
			#~ raise NotCreatedResultDir(dir_name) ;

#~ class ResultDirCreationError(Exception) :
	#~ def __init__(self,dirname="unknown") :
		




if __name__ == "__main__" :
    test_parse_gtf_content()
