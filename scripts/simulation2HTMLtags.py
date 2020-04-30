#!/usr/bin/env python3

# Internal modules
from simulation2HTMLparse import *
from simulation2HTMLplots import *



def get_html_header():
    return '''
<!doctype html>
<html lang="en">
  <head>
    <!-- Required meta tags -->
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <meta name="description" content="IntronSeeker simulation report">
    <meta name="author" content="Bardou P., Maman S., Lasguignes E., Oudin F., Cabanettes F., Klopp C.">

    <!-- Bootstrap CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css" integrity="sha384-Gn5384xqQ1aoWXA+058RXPxPg6fy4IWvTNh0E263XmFcJlSAwiGgFAW/dAiS6JXm" crossorigin="anonymous">
    <!-- Icons -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/open-iconic/1.1.1/font/css/open-iconic-bootstrap.min.css">
    <!-- Datatable -->
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.16/css/dataTables.bootstrap4.min.css">

    <!-- Personalised css -->
    <style type="text/css">
      body { font-size: .875rem; }
      .sidebar {
        position: fixed;
        top: 0;
        bottom: 0;
        left: 0;
        z-index: 100; /* Behind the navbar */
        padding: 0;
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .1);
      }
      .sidebar-sticky {
        position: -webkit-sticky;
        position: sticky;
        top: 48px; /* Height of navbar */
        height: calc(100vh - 48px);
        padding-top: 4rem;
        overflow-x: hidden;
        overflow-y: auto; /* Scrollable contents if viewport is shorter than content. */
      }
      .sidebar .nav-link {
        margin-right: 4px;
        color: #212529;
      }
      .sidebar .nav-link:hover,
      .sidebar .nav-link.active { color: #999; }
      .sidebar-heading {
        font-size: .75rem;
        text-transform: uppercase;
      }
      .navbar-brand {
        padding-top: .75rem;
        padding-bottom: .75rem;
        font-size: 1rem;
        background-color: rgba(0, 0, 0, .7);
        box-shadow: inset -1px 0 0 rgba(0, 0, 0, .7);
				z-index: 200;
      }
      .navbar .form-control {
        padding: .75rem 1rem;
        border-width: 0;
        border-radius: 0;
				z-index: 200;
      }
      .form-control-dark {
        color: #fff;
        background-color: rgba(255, 255, 255, .1);
        border-color: rgba(255, 255, 255, .1);
      }
      .form-control-dark:focus {
        border-color: transparent;
        box-shadow: 0 0 0 3px rgba(255, 255, 255, .25);
      }
      .border-top { border-top: 1px solid #e5e5e5; }
      .border-bottom { border-bottom: 1px solid #e5e5e5; }
	  .valn { vertical-align: middle !important; }
	  .anchor{
  			display: block;
  			height: 83px; /*same height as header*/
  			margin-top: -83px; /*same height as header*/
  			visibility: hidden;
	  }
    </style>

    <!-- Plotly -->
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script> 

    <title>IntronSeeker simulation report</title>
  </head>
'''

# icons : https://useiconic.com/open/
def get_html_body1(flagstat="", candidat="", assemblathon=""):
    r = '''
  <body>
    <nav class="navbar navbar-dark sticky-top bg-dark flex-md-nowrap p-0">
      <a class="navbar-brand col-sm-3 col-md-2 mr-0" href="#">IntronSeeker simulation report</a>
    </nav>

    <div class="container-fluid">
      <div class="row">
        <nav class="col-md-2 d-none d-md-block bg-light sidebar">
          <div class="sidebar-sticky">
            <ul class="nav flex-column">
              <li class="nav-item">
                <a class="nav-link" href="#inputs-parameters">
                  <span class="oi oi-file" aria-hidden="true"></span>
                  Input files
                </a>
              </li>
			  <li class="nav-item">
				<a class="nav-link" href="#ref-descr">
				  <span class="oi oi-collapse-down" aria-hidden="true"></span>
					Contigs
				</a>
			  </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#gstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Global statistics
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#nb_ctg_by_feature">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. by feature type
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#ctg_descr">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. with same features
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#contigs_len_dist">
				        <span class="oi oi-bar-chart" aria-hidden="true"></span>
				        Contigs len. distribution
				    </a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#feat_len_dist">
				        <span class="oi oi-graph" aria-hidden="true"></span>
				        Features len. distribution
				    </a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#feat_dist_intron">
				        <span class="oi oi-bar-chart" aria-hidden="true"></span>
				        Introns positions
				    </a>
			    </li>  
			  <li class="nav-item">
				<a class="nav-link" href="#read-descr">
				  <span class="oi oi-collapse-up" aria-hidden="true"></span>
					Reads
				</a>
			  </li>
			  <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#readgstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Global statistics
			    	</a>
			    </li>   
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#abundstat">
				    	<span class="oi oi-graph" aria-hidden="true"></span>
				    	Abundance
			    	</a>
			    </li>'''              
    if flagstat:
        r += '''      
              <li class="nav-item">
				<a class="nav-link" href="#flag-descr">
				  <span class="oi oi-collapse-up" aria-hidden="true"></span>
					Mapping
				</a>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#flagstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Alignment statistics
			    	</a>
			    </li>     
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#pie">
				    	<span class="oi oi-pie-chart" aria-hidden="true"></span>
				    	Mapping pie
			    	</a>
			    </li>'''      
    r += '''          
                <li class="nav-item">
				   <a class="nav-link" href="#results">
				        <span class="oi oi-folder" aria-hidden="true"></span>
					    Feature(s) extraction results
				   </a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#split">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Split detection
			    	</a>
			    </li>'''
    if assemblathon:
        r += '''                 
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#assemblystat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Intron insertion
			    	</a>
			    </li>'''      
    if candidat:
        r += '''             
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#candidatstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Candidats statistics
			    	</a>
			    </li>'''
    r +='''            
                <li class="nav-item">
                    <a class="nav-link" href="#glossary">
                        <span class="oi oi-question-mark" aria-hidden="true"></span>
                        Glossary
                    </a>
                </li>
            </ul>
          </div>
          <div style="text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px">
		    Produced by IntronSeeker_v1.0<br>
		    Copyright © 2020, <img style="width:18px;padding-bottom:2px" src="https://www.inrae.fr/themes/custom/inrae_socle/favicon.ico"><!--<img src="http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico">-->
		    <a style="color:#212529;" href="https://inrae.fr" target="_blank">INRAE</a><br>
		    Designed by the <a style="color:#212529;" href="http://sigenae.org" target="_blank">Sigenae</a> team.
          </div>
        </nav>
		<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">
'''
    return r
    

def get_html_inputfiles(files:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
			<h1 class="h4">Input files</h1>
			<span class="anchor" id="inputs-parameters"></span>
		</div>
		<div class="d-flex">
			<div class="mt-4 mr-4 pl-0 col-md-6">
				<h5>Files</h5>
				<ul class="list-group">'''
    for f in files:
        label, name = f.split("#")
        r += '''
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">'''+label+''' <span class="badge badge-success badge-pill ml-4">'''+name+'''</span></li>'''
    r += '''
				</ul>
			</div>
		</div>
'''
    return r


def get_html_seq_descr(global_stat:dict, nb_ctg_by_feature:dict, ctg_descr:dict, gtf:str, pos:dict, df_fasta:dict, df_mfasta:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Contigs</h1>
                <span class="anchor" id="ref-descr"></span>
        </div>
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Contigs statistics</h5>
                <span class="anchor" id="gstat"></span>
'''+dict_to_table(global_stat,7,True)+'''
            </div>   
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Number of sequences by feature type</h5>
                <span class="anchor" id="nb_ctg_by_feature"></span>
'''+dict_to_table(nb_ctg_by_feature,-1,False)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Number of sequences with same feature(s)</h5>
                <span class="anchor" id="nb_ctg_by_feature"></span>
'''+dict_to_table(ctg_descr,-1,False)+'''
            </div>
        </div>
'''

    # Len dist for gtf
    len_by_features, feature_names = len_dist_from_gtf(gtf)
    r += '''
        <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5>Contigs length distribution</h5>
                <span class="anchor" id="contigs_len_dist"></span>
'''+plot_hist_contigs_len(df_fasta['length'], df_mfasta['length'])+'''
        </div>
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-6">
                <h5>Features length distribution</h5>
                <span class="anchor" id="feat_len_dist"></span>
'''+plot_dist_features_len(len_by_features, feature_names)+'''
            </div>
'''
    r += '''
            <div class="mt-4 mr-0 pl-0 col-md-6">
                <h5>Insertion position of introns along contigs</h5>
                <span class="anchor" id="feat_dist_intron"></span>
'''+plot_insertion_in_contig(pos)+'''
            </div>
        </div>  
'''
    return r


def get_html_abundance(df_fasta:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5>Pourcentage of read abundance for each contig</h5>
                <span class="anchor" id="abundstat"></span>
'''+plot_abondance_model(df_fasta)+'''
            </div>
        </div>    
'''
    return r


def get_html_mapping_descr(alignments:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5> Counting table and barplots of mapped covering reads' main characteristics</h5>
                <span class="anchor" id="bam"></span>
'''+ plot_covering_reads(alignments) +'''
            </div>
        </div>    
'''
    return r


def get_html_reads_descr(global_stat_fastq : dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Reads</h1>
                <span class="anchor" id="read-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
                <h5>Global statistics</h5>
                <span class="anchor" id="readgstat"></span>
'''+dict_to_table(global_stat_fastq,-1,True)+'''
            </div>
        </div>  
'''
    return r
    
#def get_html_bam_descr(global_stat_flagstat : dict):
#    r = '''
#        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
#            <h1 class="h4">Bam (remove ?)</h1>
#                <span class="anchor" id="read-descr"></span>
#        </div>
#		<div class="d-flex">
#            <div class="mt-4 mr-0 pl-0 col-md-4">
#                <span class="anchor" id="readastat"></span>
#'''+dict_to_table(global_stat_flagstat,-1,True)+'''
#            </div>
#        </div>
#'''
#    return r    


def get_html_flagstat_descr(global_stat_flagstat:dict,df_flag_all:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Mapping</h1>
                <span class="anchor" id="flag-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Alignment statistics</h5>
                <span class="anchor" id="flagstat"></span>
'''+dict_to_table(global_stat_flagstat,-1,True)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <h5>Mapping pie</h5>
                <span class="anchor" id="pie"></span>
'''+ plot_flagstat(df_flag_all) +'''
            </div>
        </div>
'''
    return r
   
def get_html_assemblathon_descr(global_stat_assemblathon:dict):
    r = '''
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
            <h5>Assemblathon statistics</h5>
                <span class="anchor" id="assemblystat"></span>
'''+dict_to_table(global_stat_assemblathon,-1,True)+'''
            </div>
        </div>
'''
    return r

def get_html_table_descr(global_stats_table):
    r = '''
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
            <h5>TO COMPLETE</h5>
                <span class="anchor" id="assemblystat"></span>
'''+dict_to_table(global_stats_table,-1,True)+'''
            </div>
        </div>
'''
    return r


def get_html_split(mapping_bam:dict):
    r = '''
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <h5>Barplots of split reads signal detection</h5>
                <span class="anchor" id="split"></span>
''' + '''
            </div>
        </div>    
'''
    return r

    
def get_html_candidat_descr(global_stat_candidat:dict, df_candidat:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Feature(s) extraction results</h1>
                <span class="anchor" id="results"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Candidats</h5>
                <span class="anchor" id="candidatstat"></span>
'''+dict_to_table(global_stat_candidat,-1,True)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-8">
                <h5>Candidats depth distribution</h5>
                <span class="anchor" id="contigs_len_dist"></span>
'''+plot_hist_candidats_depth(df_candidat['depth'])+'''
            </div>
        </div>  
'''

    return r      
    
def get_html_glossary():
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Glossary</h1>
                <span class="anchor" id="glossary"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <span class="anchor" id="glossary"></span>
                Contig :  A contig (from contiguous) is a set of overlapping DNA segments that together represent a consensus region of DNA (source :Wikipedia)</br>
                Candidat : </br>
                Feature : </br>
                Intron : An intron is any nucleotide sequence within a gene that is removed by RNA splicing during maturation of the final RNA product. (source : Wikipedia)</br>
                Reads :  In DNA sequencing, a read is an inferred sequence of base pairs (or base pair probabilities) corresponding to all or part of a single DNA fragment.  (source : Wikipedia)</br>
                Sequence :  </br>
                Splice events : </br>
                <img src="https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/raw/master/doc/splice_event.png" alt="Splice events" style="width:825;height:245;"> </br>
                Number of features in GTF :  Number of distinct features from all GTF lines</br>
                Number of sequences by feature type :  Number of ctg by feature from all GTF lines (Ex: "Exon" see in X ctg, "Intron" see in Y ctg, ...)</br>
                Number of sequences with same feature(s) : Number of features profiles by ctg (Ex: "1 Exon & 2 Intron" see in X ctg, "3 Introns" see in Y ctg , ...) </br>
            </div>
        </div>  
'''

    return r    

def get_html_footer():
    return '''
        </main>
      </div>
    </div>

    <!-- jQuery first, then Popper.js, then Bootstrap JS -->
    <script type="text/javascript" language="javascript" src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
    <script type="text/javascript" language="javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

    <!-- Datatable -->
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js"></script>
  </body>
</html>
'''


# Return HTML table from dict
# Param 1 : dict (key, val)
# Param 2 : int for which line number first col (key) will be right align
# Param 3 : bool to remove first char of the key (used to sort by key)
def dict_to_table(d : dict, i : int, rmfirstchar : bool):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    c = 0
    for k, v in sorted(d.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c >= i & i!=-1):
            table += " text-right"
        val=re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str(v))   # Remove (nb%) in val    
        if(rmfirstchar):
            table += "'>" + k[1:] + "</td><td class='valn text-right'>" + split_int(round(float(val)), ' ') + "</td></tr>"
        else:
            table += "'>" + k + "</td><td class='valn text-right'>" + split_int(round(float(val)), ' ') + "</td></tr>"
        c += 1
    table += "</tbody></table>"
    return table

# Return HTML table from dict
# Param 1 : dict (key, val)
# Param 2 : int for which line number first col (key) will be right align
# Param 3 : bool to remove first char of the key (used to sort by key)
def dict2_to_table(d1 : dict, d2:dict, i : int, rmfirstchar : bool):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    c = 0
    c2 = 0
    for k, v in sorted(d1.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c >= i & i!=-1):
            table += " text-right"
        if(rmfirstchar):
            table += "'>" + k[1:] + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        else:
            table += "'>" + k + "</td><td class='valn text-right'>" + str(v) + "</td></tr>"
        c += 1
    for k2, v2 in sorted(d2.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn"
        if(c2 >= i & i!=-1):
            table += " text-right"
        if(rmfirstchar):
            table += "'>" + k2[1:] + "</td><td class='valn text-right'>" + str(v2) + "</td></tr>"
        else:
            table += "'>" + k2 + "</td><td class='valn text-right'>" + str(v2) + "</td></tr>"
        c2 += 1    
    table += "</tbody></table>"
    return table