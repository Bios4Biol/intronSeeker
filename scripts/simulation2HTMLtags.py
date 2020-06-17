#!/usr/bin/env python3

# Internal modules
from simulation2HTMLparse import *
from simulation2HTMLplots import *
from collections import defaultdict


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
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/open-iconic/1.1.1/font/css/open-iconic-bootstrap.min.css ">
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
def get_html_body1(flagstat="", split="", candidat=""):
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
				    	Contigs statistics
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#ref-descr">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. by feature type
			    	</a>
			    </li>
			    <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#ref-descr">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Nb. of seq. with same features
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#ref-descr_assemblathon_fasta">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Assemblathon(s)
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
				    	Reads statistics
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
			   </li>'''
    r += '''          
                <li class="nav-item">
				   <a class="nav-link" href="#features-descr">
				        <span class="oi oi-folder" aria-hidden="true"></span>
					    Intron extraction results
				   </a>
			    </li>'''
    if split:
        r += '''    
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#splitstat">
				    	<span class="oi oi-graph" aria-hidden="true"></span>
				    	Split reads
			    	</a>
			    </li>'''                         
    if candidat:
        r += '''             
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#candidatstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Detected introns
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#filtered_detected_introns">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Filtered detected introns
			    	</a>
			    </li>
                '''
    r +='''      
                <li class="nav-item">
                    <a class="nav-link" href="#precision">
                        <span class="oi oi-circle-check" aria-hidden="true"></span>
                        Comparison features / detected introns
                    </a>
                </li> 
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#precision">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	GTF vs detected introns
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#precision">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Filtered features
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#precision">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Confusion matrix
			    	</a>
			    </li>
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#precision">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Precision, recall & F1 score
			    	</a>
			    </li>
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
		    Copyright © 2020, <img style="width:18px;padding-bottom:2px" src="https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/raw/master/scripts/icon/favicon.ico"   ><!--<img src="http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico">-->
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


def get_html_seq_descr(global_stat:dict, nb_ctg_by_feature:dict, ctg_descr:dict, gtf:str, pos:dict, df_fasta:dict, df_mfasta:dict, global_stat_assemblathon_fasta:dict, global_stat_assemblathon_mfasta:" "):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Contigs</h1>
                <span class="anchor" id="ref-descr"></span>
        </div>
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <div class="d-flex">
                    <div>
                        <h5>Contigs statistics</h5>
                        <span class="anchor" id="ref-descr"></span>
                    </div>
                    <div class="ml-2">
                        <span
                            class="badge badge-pill badge-dark"
                            data-toggle="tooltip"
                            data-placement="top"
                            data-html="true"
                            title="<p nowrap><b>Number of features in GTF</b>: total number of features from GTF file</p><p nowrap><b>Number of distinct features in GTF</b>: number of distinct features from all GTF lines</p>">
                            <span class="oi oi-info" ></span>
                        </span>
                    </div>
                </div>               
'''+dict_to_table(global_stat,7,True)+'''
            </div>   
            <div class="mt-4 mr-0 pl-0 col-md-4">
                <div class="d-flex">
                    <div>
                        <h5>Number of sequences by feature type</h5>
                        <span class="anchor" id="ref-descr"></span>
                    </div>
                    <div class="ml-2">
                        <span
                            class="badge badge-pill badge-dark"
                            data-toggle="tooltip"
                            data-placement="top"
                            data-html="true"
                            title="<p nowrap><b>Number of ctg by feature from all GTF lines </b> (Ex: 'Exon' see in X ctg, 'Intron' see in Y ctg, ...)</p>">
                            <span class="oi oi-info"></span>
                        </span>
                    </div>
                </div>
'''+dict_to_table(nb_ctg_by_feature,-1,False)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
               <div class="d-flex">
                    <div>
                        <h5>Number of sequences with same feature(s)</h5>
                        <span class="anchor" id="ref-descr"></span>
                    </div>
                    <div class="ml-2">
                        <span
                            class="badge badge-pill badge-dark"
                            data-toggle="tooltip"
                            data-placement="top"
                            data-html="true"
                            title="<p nowrap><b>Number of features profiles by ctg </b>(Ex: '1 Exon & 2 Intron' see in X ctg, '3 Introns' see in Y ctg , ...)</p>">
                            <span class="oi oi-info"></span>
                        </span>
                    </div>
                </div>
'''+dict_to_table(ctg_descr,-1,False)+'''
            </div>
        </div>
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-6">
                <h5>Assemblathon statistics</h5>
                <span class="anchor" id="ref-descr_assemblathon_fasta"></span>
'''+dict_to_table_multi_col("Contig FASTA", global_stat_assemblathon_fasta, "Contig FASTA with feature(s)", global_stat_assemblathon_mfasta)+'''
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


def get_html_reads_descr(global_stat_fastq : dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Reads</h1>
                <span class="anchor" id="read-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
                <h5>Reads statistics</h5>
                <span class="anchor" id="readgstat"></span>
'''+dict_to_table(global_stat_fastq,-1,True)+'''
            </div>
        </div>  
'''
    return r
    
def get_html_flagstat_descr(global_stat_flagstat:dict):
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Mapping</h1>
                <span class="anchor" id="flag-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
'''+ dict_to_table(global_stat_flagstat, -1,True) +'''
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


# Results
def get_html_results():
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Intron extraction results</h1>
                <span class="anchor" id="features-descr"></span>
        </div>
'''
    return r   


# Plots split detection stats   
def get_html_split_descr(df_splitRead:dict):
    r =  '''
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Split reads</h5>
                <span class="anchor" id="splitstat"></span>
''' +  dict_to_table(df_splitRead, 2, True) + '''  
            </div>
        </div>  
'''
    return r   


# Plots detected candidats stats : tab and graph
def get_html_detected(global_stat_detected_introns:dict, df_candidat:dict):
    r =  '''
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Detected introns</h5>
                <span class="anchor" id="candidatstat"></span>
'''+dict_to_table(global_stat_detected_introns,-1,True)+'''
            </div>  
            <div class="mt-4 mr-0 pl-0 col-md-8">
                <h5>Detected introns depth distribution</h5>
                <span class="anchor" id="candidatstat"></span>
'''+plot_hist_candidats_depth(df_candidat['depth'])+'''
            </div>
        </div>  

'''
    return r   


# Tab candidat filtered detected candidats (and if simulation detectable features)
def get_html_candidat(global_stat_f_detected_introns:dict, global_stat_detectable_features=0):
    r =  '''
        <div class="d-flex"> 
            <div class="mt-4 mr-0 pl-0'''
    if global_stat_detectable_features:
        r += ''' col-md-8">
            <h5>Filtered detected introns and filtered features (detectable)</h5>
                <span class="anchor" id="filtered_detected_introns"></span>'''
        r += dict_to_table_multi_col("Filtered detected introns", global_stat_f_detected_introns,
                                     "Filtered features (detectable)", global_stat_detectable_features)
    else:
        r += ''' col-md-4">
            <h5>Filtered detected introns</h5>
                <span class="anchor" id="filtered_detected_introns"></span>'''
        r += dict_to_table(global_stat_f_detected_introns, -1, True)
    r+='''
            </div>
        </div>
'''
    return r   


# precision, recall and F1 score   
def get_html_precision(global_stat_precision:dict, TP:int, TN:int, FP:int, FN:int):
    r =  '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Comparison between features and detected introns, specificity and sensibility</h1>
                <span class="anchor" id="precision"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Confusion matrix</h5>
                <span class="anchor" id="precision"></span>
'''+dict_to_table_3col(TP, TN, FP, FN)+'''
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Precision, recall and F1 score</h5>
                <span class="anchor" id="precision"></span>
'''+dict_to_table(global_stat_precision,-1,True)+'''
            </div>
        </div>  
'''
    return r       


# Define all term used in this simulation report
def get_html_glossary():
    r = '''
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Glossary</h1>
                <span class="anchor" id="glossary"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <span class="anchor" id="glossary"></span>
                <br>
                <img src="https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/raw/master/doc/IntronSeeker-glossary.png" alt="intronSeeker glossary" style="width:825;height:245;">
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
    <script type="text/javascript" language="javascript" src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js " integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>

    <!-- Datatable -->
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script type="text/javascript" language="javascript" src="https://cdn.datatables.net/1.10.16/js/dataTables.bootstrap4.min.js"></script>
    <script type="text/javascript" language="javascript">
        $(function () {
            $('[data-toggle="tooltip"]').tooltip()
        })
    </script>
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
        val = re.sub(r'(\([a-zA-Z0-9_]*.[a-zA-Z0-9_]*%\))', r" ", str(v))   # Remove (nb%) in val    
        values = split_int(round(float(val)), ' ')   
        if(rmfirstchar):
            table += "'>" + k[1:] + "</td><td class='valn text-right'>" + values + "</td></tr>"   
        else:
            table += "'>" + k + "</td><td class='valn text-right'>" + values + "</td></tr>"
        c += 1
    table += "</tbody></table>"
    return table


def dict_to_table_multi_col(col1name:str, d1:dict, col2name="", d2="", col3name="", d3=""):
    colmd = "width='50%'"
    if col2name:
        colmd = "width='33%'"
    if col3name:
        colmd = "width='25%'"
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <thead>
''' 
    table += "<tr><td " + colmd + "></td><td class='valn text-center' " + colmd + ">" + col1name + "</td>"
    if col2name:
        table += "<td class='valn text-center' " + colmd + ">" + col2name + "</td>"
    if col3name:
        table += "<td class='valn text-center' " + colmd + ">" + col3name + "</td>"
    table += "</tr></thead><tbody>"
    for k, v in sorted(d1.items(), key=lambda t: t[0]):
        table += "<tr><td class='valn text-left'>" + k[1:] + "</td>"
        table += "<td class='valn text-right'>" + split_int(round(float(v)), ' ')  + "</td>"
        if col2name:
            table += "<td class='valn text-right'>" + split_int(round(float(d2[k])), ' ')  + "</td>"
        if col3name:
            table += "<td class='valn text-right'>" + split_int(round(float(d3[k])), ' ')  + "</td>"
        table += "</tr>"
    return table + "</tbody></table>"    


# Return HTML table from dict (as dict_to_table function) but with 3 columns instead of 2
# Param 1 : TP is the number of detectable and found features (int value)
# Param 2 : TN is the number of detectable and not found features (int value)
# Param 3 : FP is the number of undetectable and found features (int value)
# Param 4 : FN is the number of undetectable and not found features (int value)
def dict_to_table_3col(TP: int, TN: int, FP: int, FN:int):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    table += "<tr> <td></td><td class='valn text-center'>Not found</td><td class='valn text-center'>Found</td></tr>"
    table += "<tr><td class='valn text-center'>Undetectable</td><td class='valn text-center'>" + str(FN) + "</td><td class='valn text-center'>" + str(FP) + "</td><td class='valn text-right'>Actual detectable</td></tr>"  
    table += "<tr><td class='valn text-center'>Detectable</td><td class='valn text-center'>" + str(TN) + "</td><td class='valn text-center'>" + str(TP) + "</td><td class='valn text-right'>"+str(TN+TP)+"</td></tr>"  
    table += "<tr><td></td><td class='valn text-right'>Total detectable found</td><td class='valn text-center'>" + str(TP+FP) + "</td></tr>"  
    table += "</tbody></table>"
    return table

# Return HTML table from dataframe
# Param 1 : dict (key, val)
# Param 2 : int for which line number first col (key) will be right align
# Param 3 : bool to change text side in cells' tab
def df_to_table(df : dict, i : int, changeTextSide : bool):
    table = '''
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
'''
    j = 0
    for row in df.itertuples():
        if(j >= i & changeTextSide):
            table += "<tr><td class='valn text-right'>" + str(row.titles) + "</td><td class='valn text-right'>" + str(split_int(round(row.values))) + "</td></tr>"
        else:
            table += "<tr><td class='valn text-left'>" + str(row.titles) + "</td><td class='valn text-right'>" + str(split_int(round(row.values))) + "</td></tr>"
        j += 1
    
    table += "</tbody></table>"
    return table