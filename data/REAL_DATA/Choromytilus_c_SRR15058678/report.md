
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
				    <a class="nav-link" href="#ref-descr">
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
              <li class="nav-item">
				<a class="nav-link" href="#flag-descr">
				  <span class="oi oi-collapse-up" aria-hidden="true"></span>
					Mapping
				</a>
			   </li>          
                <li class="nav-item">
				   <a class="nav-link" href="#features-descr">
				        <span class="oi oi-folder" aria-hidden="true"></span>
					    Intron extraction results
				   </a>
			    </li>    
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#splitstat">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Split reads
			    	</a>
			    </li>             
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
                <li class="nav-item" style="padding-left:10px">
				    <a class="nav-link" href="#too_complex_detected">
				    	<span class="oi oi-list" aria-hidden="true"></span>
				    	Contigs with too complex introns
			    	</a>
			    </li>      
                <li class="nav-item">
                    <a class="nav-link" href="#eval">
                        <span class="oi oi-circle-check" aria-hidden="true"></span>
                        Evaluation of the detected introns (simulation only)
                    </a>
                </li> 
                <li class="nav-item">
                    <a class="nav-link" href="#glossary">
                        <span class="oi oi-info" aria-hidden="true"></span>
                        Glossary
                    </a>
                </li>
            </ul>
          </div>
          <div style="text-align:center;font-size:smaller;color:darkgrey;margin-top:-25px">
		    Produced by IntronSeeker_v1.0<br>
		    Copyright © 2020, <img style="width:18px;padding-bottom:2px" src="https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/raw/master/scripts/icon/favicon.ico"   ><!--<img src="http://www.inra.fr/extension/itkinra/design/inra/images/favicon.ico">-->
		    <a style="color:#212529;" href="https://www.inrae.fr" target="_blank">INRAE</a><br>
		    Designed by the <a style="color:#212529;" href="http://sigenae.org" target="_blank">Sigenae</a> team.
          </div>
        </nav>
		<main role="main" class="col-md-9 ml-sm-auto col-lg-10 pt-3 px-4">

        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pb-2 border-bottom">
			<h1 class="h4">Input files</h1>
			<span class="anchor" id="inputs-parameters"></span>
		</div>
		<div class="d-flex">
			<div class="mt-4 mr-4 pl-0 col-md-6">
				<h5>Files</h5>
				<ul class="list-group">
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Contig FASTA <span class="badge badge-success badge-pill ml-4">GJJD01.1.fsa_nt</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Read1 FASTQ <span class="badge badge-success badge-pill ml-4">SRR15058678_1.fastq</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Flagstat <span class="badge badge-success badge-pill ml-4">hisat2_GJJD01.sort.flagstat.txt</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Split <span class="badge badge-success badge-pill ml-4">srs_GJJD01_split_alignments.txt</span></li>
                    <li class="list-group-item d-flex justify-content-between align-items-center p-2">Candidat <span class="badge badge-success badge-pill ml-4">srs_GJJD01_candidates.txt</span></li>
				</ul>
			</div>
		</div>

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

            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Contig FASTA - Number of sequences</td><td class='valn text-right'>98 794</td></tr><tr><td class='valn'>Contig FASTA - Mean sequence length</td><td class='valn text-right'>827</td></tr></tbody></table>	
                </div>   
        </div>

        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Reads</h1>
                <span class="anchor" id="read-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-4 pl-0 col-md-4">
                 <h5>Reads statistics</h5>
                <span class="anchor" id="readgstat"></span>
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Number of reads</td><td class='valn text-right'>10 490 833</td></tr><tr><td class='valn'>Mean coverage</td><td class='valn text-right'>12.92</td></tr><tr><td class='valn'>Min reads length</td><td class='valn text-right'>15</td></tr><tr><td class='valn'>Max reads length</td><td class='valn text-right'>101</td></tr><tr><td class='valn'>Mean reads length</td><td class='valn text-right'>100.59</td></tr></tbody></table>
            </div>
            <div class="mt-4 mr-0 pl-0 col-md-8">
                    <h5>Reads length distribution</h5>
                    <span class="anchor" id="readgstat"></span><div>                            <div id="22049aa2-3c94-4b63-ad05-b871beca287b" class="plotly-graph-div" style="height:100%; width:100%;"></div>            <script type="text/javascript">                                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("22049aa2-3c94-4b63-ad05-b871beca287b")) {                    Plotly.newPlot(                        "22049aa2-3c94-4b63-ad05-b871beca287b",                        [{"name":"Reads length distribution","opacity":0.85,"x":[15.0,15.172,15.344,15.516,15.688,15.86,16.032,16.204,16.376,16.548,16.72,16.892,17.064,17.236,17.408,17.58,17.752,17.924,18.096,18.268,18.439999999999998,18.612,18.784,18.956,19.128,19.3,19.472,19.644,19.816,19.988,20.16,20.332,20.503999999999998,20.676,20.848,21.02,21.192,21.364,21.536,21.708,21.88,22.052,22.224,22.396,22.567999999999998,22.74,22.912,23.084,23.256,23.427999999999997,23.6,23.772,23.944,24.116,24.287999999999997,24.46,24.631999999999998,24.804,24.976,25.148,25.32,25.491999999999997,25.664,25.836,26.008,26.18,26.351999999999997,26.524,26.695999999999998,26.868,27.04,27.212,27.384,27.555999999999997,27.728,27.9,28.072,28.244,28.415999999999997,28.588,28.759999999999998,28.932,29.104,29.275999999999996,29.448,29.619999999999997,29.791999999999998,29.964,30.136,30.308,30.479999999999997,30.652,30.823999999999998,30.996,31.168,31.34,31.512,31.683999999999997,31.855999999999998,32.028,32.2,32.372,32.544,32.715999999999994,32.888,33.06,33.232,33.403999999999996,33.57599999999999,33.748,33.92,34.092,34.263999999999996,34.436,34.608,34.78,34.952,35.123999999999995,35.296,35.468,35.64,35.812,35.983999999999995,36.156,36.328,36.5,36.672,36.843999999999994,37.016,37.188,37.36,37.532,37.70399999999999,37.876,38.048,38.22,38.391999999999996,38.56399999999999,38.736,38.908,39.08,39.251999999999995,39.424,39.596,39.768,39.94,40.111999999999995,40.284,40.456,40.628,40.8,40.971999999999994,41.144,41.316,41.488,41.66,41.831999999999994,42.004,42.176,42.348,42.519999999999996,42.69199999999999,42.864,43.036,43.208,43.379999999999995,43.55199999999999,43.724,43.896,44.068,44.239999999999995,44.412,44.583999999999996,44.756,44.928,45.099999999999994,45.272,45.444,45.616,45.788,45.959999999999994,46.132,46.304,46.476,46.647999999999996,46.81999999999999,46.992,47.163999999999994,47.336,47.507999999999996,47.68,47.852,48.024,48.196,48.367999999999995,48.54,48.711999999999996,48.884,49.056,49.227999999999994,49.4,49.571999999999996,49.744,49.916,50.087999999999994,50.26,50.431999999999995,50.604,50.775999999999996,50.948,51.12,51.291999999999994,51.464,51.635999999999996,51.808,51.98,52.151999999999994,52.324,52.495999999999995,52.668,52.839999999999996,53.012,53.184,53.355999999999995,53.528,53.699999999999996,53.872,54.044,54.215999999999994,54.388,54.559999999999995,54.732,54.903999999999996,55.07599999999999,55.248,55.419999999999995,55.592,55.763999999999996,55.936,56.108,56.279999999999994,56.452,56.623999999999995,56.796,56.967999999999996,57.13999999999999,57.312,57.483999999999995,57.656,57.827999999999996,58.0,58.172,58.343999999999994,58.516,58.687999999999995,58.86,59.032,59.20399999999999,59.376,59.547999999999995,59.72,59.891999999999996,60.06399999999999,60.236,60.407999999999994,60.58,60.751999999999995,60.924,61.096,61.267999999999994,61.44,61.611999999999995,61.784,61.955999999999996,62.12799999999999,62.3,62.471999999999994,62.644,62.815999999999995,62.988,63.16,63.331999999999994,63.504,63.675999999999995,63.848,64.02,64.192,64.364,64.536,64.708,64.88,65.05199999999999,65.22399999999999,65.39599999999999,65.568,65.74,65.912,66.084,66.256,66.428,66.6,66.77199999999999,66.94399999999999,67.11599999999999,67.288,67.46,67.632,67.804,67.976,68.148,68.32,68.49199999999999,68.66399999999999,68.836,69.008,69.17999999999999,69.352,69.524,69.696,69.868,70.03999999999999,70.21199999999999,70.38399999999999,70.556,70.728,70.9,71.072,71.244,71.416,71.588,71.75999999999999,71.93199999999999,72.10399999999998,72.276,72.448,72.62,72.792,72.964,73.136,73.30799999999999,73.47999999999999,73.65199999999999,73.824,73.996,74.16799999999999,74.34,74.512,74.684,74.856,75.02799999999999,75.19999999999999,75.37199999999999,75.544,75.716,75.888,76.06,76.232,76.404,76.576,76.74799999999999,76.91999999999999,77.09199999999998,77.264,77.43599999999999,77.608,77.78,77.952,78.124,78.29599999999999,78.46799999999999,78.63999999999999,78.812,78.984,79.15599999999999,79.32799999999999,79.5,79.672,79.844,80.01599999999999,80.18799999999999,80.36,80.532,80.704,80.87599999999999,81.048,81.22,81.392,81.564,81.73599999999999,81.908,82.08,82.252,82.42399999999999,82.59599999999999,82.768,82.94,83.112,83.28399999999999,83.45599999999999,83.628,83.8,83.972,84.14399999999999,84.31599999999999,84.488,84.66,84.832,85.00399999999999,85.17599999999999,85.348,85.52,85.692,85.86399999999999,86.036,86.208,86.38,86.55199999999999,86.72399999999999,86.896,87.068,87.24,87.41199999999999,87.58399999999999,87.756,87.928,88.1,88.27199999999999,88.44399999999999,88.616,88.788,88.96,89.13199999999999,89.30399999999999,89.476,89.648,89.82,89.99199999999999,90.16399999999999,90.336,90.508,90.67999999999999,90.85199999999999,91.024,91.196,91.368,91.53999999999999,91.71199999999999,91.884,92.056,92.228,92.39999999999999,92.57199999999999,92.744,92.916,93.088,93.25999999999999,93.43199999999999,93.604,93.776,93.948,94.11999999999999,94.29199999999999,94.464,94.636,94.80799999999999,94.97999999999999,95.15199999999999,95.324,95.496,95.66799999999999,95.83999999999999,96.012,96.184,96.356,96.52799999999999,96.69999999999999,96.872,97.044,97.216,97.38799999999999,97.55999999999999,97.732,97.904,98.076,98.24799999999999,98.41999999999999,98.592,98.764,98.93599999999999,99.10799999999999,99.27999999999999,99.452,99.624,99.79599999999999,99.96799999999999,100.13999999999999,100.312,100.484,100.65599999999999,100.82799999999999],"y":[149,0,0,0,0,160,0,0,0,0,0,133,0,0,0,0,0,122,0,0,0,0,0,139,0,0,0,0,0,112,0,0,0,0,94,0,0,0,0,0,114,0,0,0,0,0,100,0,0,0,0,0,101,0,0,0,0,0,105,0,0,0,0,100,0,0,0,0,0,85,0,0,0,0,0,102,0,0,0,0,0,87,0,0,0,0,0,103,0,0,0,0,0,94,0,0,0,0,77,0,0,0,0,0,83,0,0,0,0,0,93,0,0,0,0,0,87,0,0,0,0,0,159,0,0,0,0,103,0,0,0,0,0,77,0,0,0,0,0,85,0,0,0,0,0,92,0,0,0,0,0,234,0,0,0,0,98,0,0,0,0,0,84,0,0,0,0,0,91,0,0,0,0,0,109,0,0,0,0,0,204,0,0,0,0,0,92,0,0,0,0,229,0,0,0,0,0,113,0,0,0,0,0,94,0,0,0,0,0,92,0,0,0,0,0,93,0,0,0,0,90,0,0,0,0,0,139,0,0,0,0,0,127,0,0,0,0,0,124,0,0,0,0,0,133,0,0,0,0,0,161,0,0,0,0,158,0,0,0,0,0,129,0,0,0,0,0,133,0,0,0,0,0,126,0,0,0,0,0,180,0,0,0,0,216,0,0,0,0,0,244,0,0,0,0,0,254,0,0,0,0,0,237,0,0,0,0,0,257,0,0,0,0,387,0,0,0,0,0,436,0,0,0,0,0,445,0,0,0,0,0,459,0,0,0,0,0,515,0,0,0,0,0,520,0,0,0,0,684,0,0,0,0,0,608,0,0,0,0,0,488,0,0,0,0,0,637,0,0,0,0,0,776,0,0,0,0,725,0,0,0,0,0,1041,0,0,0,0,0,1547,0,0,0,0,0,1232,0,0,0,0,0,2224,0,0,0,0,3357,0,0,0,0,0,3200,0,0,0,0,0,7030,0,0,0,0,0,9184,0,0,0,0,0,7458,0,0,0,0,0,17558,0,0,0,0,25646,0,0,0,0,0,16209,0,0,0,0,0,36416,0,0,0,0,0,60851,0,0,0,0,0,17404,0,0,0,0,54475,0,0,0,0,0,236995,0,0,0,0,0,51091,0,0,0,0,0,64815,0,0,0,0,0,337189,0,0,0,0,9522234],"type":"bar"}],                        {"margin":{"b":20,"l":50,"pad":0,"r":50,"t":30},"template":{"data":{"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"#E5ECF6","showlakes":true,"showland":true,"subunitcolor":"white"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2}}},"xaxis":{"title":{"text":"Length (bp)"}},"yaxis":{"title":{"text":"Number of reads"}}},                        {"responsive": true}                    )                };                            </script>        </div>
            </div>
        </div>
        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Mapping (QC-passed reads)</h1>
                <span class="anchor" id="flag-descr"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">

            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Number of reads</td><td class='valn text-right'>0</td></tr><tr><td class='valn'>Number of lines in the BAM</td><td class='valn text-right'>10 526 093</td></tr><tr><td class='valn'>Number of mapped reads</td><td class='valn text-right'>7 487 211</td></tr><tr><td class='valn'>Percentage of mapped reads</td><td class='valn text-right'>71.37</td></tr><tr><td class='valn'>Number of properly paired reads</td><td class='valn text-right'>0</td></tr><tr><td class='valn'>Percentage of properly paired reads</td><td class='valn text-right'>0</td></tr><tr><td class='valn'>Secondary alignments</td><td class='valn text-right'>35 260</td></tr><tr><td class='valn'>Singletons alignements</td><td class='valn text-right'>0</td></tr></tbody></table>
            </div>
        </div>

        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Intron extraction results</h1>
                <span class="anchor" id="features-descr"></span>
        </div>

		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Split reads</h5>
                <span class="anchor" id="splitstat"></span>

            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Number of reads overlapping introns</td><td class='valn text-right'>16 729</td></tr><tr><td class='valn'>Mean length of introns</td><td class='valn text-right'>270.77</td></tr><tr><td class='valn text-right'>Junction TC_GG</td><td class='valn text-right'>927</td></tr><tr><td class='valn text-right'>Junction NN_NN</td><td class='valn text-right'>699</td></tr><tr><td class='valn text-right'>Junction CT_GC</td><td class='valn text-right'>403</td></tr><tr><td class='valn text-right'>Junction AT_AC</td><td class='valn text-right'>338</td></tr><tr><td class='valn text-right'>Junction CT_NN</td><td class='valn text-right'>242</td></tr><tr><td class='valn text-right'>Other junctions</td><td class='valn text-right'>4 449</td></tr><tr><td class='valn text-right'>Canonical junction (GT_AG or CT_AC)</td><td class='valn text-right'>9 671</td></tr></tbody></table>  
            </div>
        </div>  

		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Detected introns</h5>
                <span class="anchor" id="candidatstat"></span>

            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Number</td><td class='valn text-right'>3 189</td></tr><tr><td class='valn'>Min length</td><td class='valn text-right'>20</td></tr><tr><td class='valn'>Max length</td><td class='valn text-right'>5 511</td></tr><tr><td class='valn'>Mean length</td><td class='valn text-right'>223.13</td></tr><tr><td class='valn'>Min depth</td><td class='valn text-right'>1</td></tr><tr><td class='valn'>Max depth</td><td class='valn text-right'>2 376</td></tr><tr><td class='valn'>Mean depth</td><td class='valn text-right'>5.25</td></tr></tbody></table>
            </div>  
            <div class="mt-4 mr-0 pl-0 col-md-8">
                <h5>Detected introns depth distribution</h5>
                <span class="anchor" id="candidatstat"></span>
<div>                            <div id="6a27e32c-8054-47f0-9f8c-eab948f23009" class="plotly-graph-div" style="height:100%; width:100%;"></div>            <script type="text/javascript">                                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("6a27e32c-8054-47f0-9f8c-eab948f23009")) {                    Plotly.newPlot(                        "6a27e32c-8054-47f0-9f8c-eab948f23009",                        [{"name":"Candidats depth","opacity":0.85,"x":[1.0,5.75,10.5,15.25,20.0,24.75,29.5,34.25,39.0,43.75,48.5,53.25,58.0,62.75,67.5,72.25,77.0,81.75,86.5,91.25,96.0,100.75,105.5,110.25,115.0,119.75,124.5,129.25,134.0,138.75,143.5,148.25,153.0,157.75,162.5,167.25,172.0,176.75,181.5,186.25,191.0,195.75,200.5,205.25,210.0,214.75,219.5,224.25,229.0,233.75,238.5,243.25,248.0,252.75,257.5,262.25,267.0,271.75,276.5,281.25,286.0,290.75,295.5,300.25,305.0,309.75,314.5,319.25,324.0,328.75,333.5,338.25,343.0,347.75,352.5,357.25,362.0,366.75,371.5,376.25,381.0,385.75,390.5,395.25,400.0,404.75,409.5,414.25,419.0,423.75,428.5,433.25,438.0,442.75,447.5,452.25,457.0,461.75,466.5,471.25,476.0,480.75,485.5,490.25,495.0,499.75,504.5,509.25,514.0,518.75,523.5,528.25,533.0,537.75,542.5,547.25,552.0,556.75,561.5,566.25,571.0,575.75,580.5,585.25,590.0,594.75,599.5,604.25,609.0,613.75,618.5,623.25,628.0,632.75,637.5,642.25,647.0,651.75,656.5,661.25,666.0,670.75,675.5,680.25,685.0,689.75,694.5,699.25,704.0,708.75,713.5,718.25,723.0,727.75,732.5,737.25,742.0,746.75,751.5,756.25,761.0,765.75,770.5,775.25,780.0,784.75,789.5,794.25,799.0,803.75,808.5,813.25,818.0,822.75,827.5,832.25,837.0,841.75,846.5,851.25,856.0,860.75,865.5,870.25,875.0,879.75,884.5,889.25,894.0,898.75,903.5,908.25,913.0,917.75,922.5,927.25,932.0,936.75,941.5,946.25,951.0,955.75,960.5,965.25,970.0,974.75,979.5,984.25,989.0,993.75,998.5,1003.25,1008.0,1012.75,1017.5,1022.25,1027.0,1031.75,1036.5,1041.25,1046.0,1050.75,1055.5,1060.25,1065.0,1069.75,1074.5,1079.25,1084.0,1088.75,1093.5,1098.25,1103.0,1107.75,1112.5,1117.25,1122.0,1126.75,1131.5,1136.25,1141.0,1145.75,1150.5,1155.25,1160.0,1164.75,1169.5,1174.25,1179.0,1183.75,1188.5,1193.25,1198.0,1202.75,1207.5,1212.25,1217.0,1221.75,1226.5,1231.25,1236.0,1240.75,1245.5,1250.25,1255.0,1259.75,1264.5,1269.25,1274.0,1278.75,1283.5,1288.25,1293.0,1297.75,1302.5,1307.25,1312.0,1316.75,1321.5,1326.25,1331.0,1335.75,1340.5,1345.25,1350.0,1354.75,1359.5,1364.25,1369.0,1373.75,1378.5,1383.25,1388.0,1392.75,1397.5,1402.25,1407.0,1411.75,1416.5,1421.25,1426.0,1430.75,1435.5,1440.25,1445.0,1449.75,1454.5,1459.25,1464.0,1468.75,1473.5,1478.25,1483.0,1487.75,1492.5,1497.25,1502.0,1506.75,1511.5,1516.25,1521.0,1525.75,1530.5,1535.25,1540.0,1544.75,1549.5,1554.25,1559.0,1563.75,1568.5,1573.25,1578.0,1582.75,1587.5,1592.25,1597.0,1601.75,1606.5,1611.25,1616.0,1620.75,1625.5,1630.25,1635.0,1639.75,1644.5,1649.25,1654.0,1658.75,1663.5,1668.25,1673.0,1677.75,1682.5,1687.25,1692.0,1696.75,1701.5,1706.25,1711.0,1715.75,1720.5,1725.25,1730.0,1734.75,1739.5,1744.25,1749.0,1753.75,1758.5,1763.25,1768.0,1772.75,1777.5,1782.25,1787.0,1791.75,1796.5,1801.25,1806.0,1810.75,1815.5,1820.25,1825.0,1829.75,1834.5,1839.25,1844.0,1848.75,1853.5,1858.25,1863.0,1867.75,1872.5,1877.25,1882.0,1886.75,1891.5,1896.25,1901.0,1905.75,1910.5,1915.25,1920.0,1924.75,1929.5,1934.25,1939.0,1943.75,1948.5,1953.25,1958.0,1962.75,1967.5,1972.25,1977.0,1981.75,1986.5,1991.25,1996.0,2000.75,2005.5,2010.25,2015.0,2019.75,2024.5,2029.25,2034.0,2038.75,2043.5,2048.25,2053.0,2057.75,2062.5,2067.25,2072.0,2076.75,2081.5,2086.25,2091.0,2095.75,2100.5,2105.25,2110.0,2114.75,2119.5,2124.25,2129.0,2133.75,2138.5,2143.25,2148.0,2152.75,2157.5,2162.25,2167.0,2171.75,2176.5,2181.25,2186.0,2190.75,2195.5,2200.25,2205.0,2209.75,2214.5,2219.25,2224.0,2228.75,2233.5,2238.25,2243.0,2247.75,2252.5,2257.25,2262.0,2266.75,2271.5,2276.25,2281.0,2285.75,2290.5,2295.25,2300.0,2304.75,2309.5,2314.25,2319.0,2323.75,2328.5,2333.25,2338.0,2342.75,2347.5,2352.25,2357.0,2361.75,2366.5,2371.25],"y":[2834,162,71,35,14,15,6,7,8,6,4,4,2,0,0,2,2,0,0,2,0,0,1,0,1,1,2,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],"type":"bar"}],                        {"margin":{"b":20,"l":50,"pad":0,"r":50,"t":30},"template":{"data":{"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"#E5ECF6","showlakes":true,"showland":true,"subunitcolor":"white"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","gridwidth":2,"linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","zerolinewidth":2}}},"xaxis":{"title":{"text":"Detected introns depth"}},"yaxis":{"title":{"text":"Number of detected introns"}}},                        {"responsive": true}                    )                };                            </script>        </div>
            </div>
        </div>  


        <div class="d-flex"> 
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Filtered detected introns</h5>
                <span class="anchor" id="filtered_detected_introns"></span>
            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>Number</td><td class='valn text-right'>109</td></tr><tr><td class='valn'>Filtered because of depth (<= 1)</td><td class='valn text-right'>2 775</td></tr><tr><td class='valn'>Filtered because of length (>= 80%)</td><td class='valn text-right'>52</td></tr><tr><td class='valn'>Filtered because of non canonical junction</td><td class='valn text-right'>1 887</td></tr><tr><td class='valn'>Filtered because of overlapping introns</td><td class='valn text-right'>0</td></tr><tr><td class='valn'>Min length</td><td class='valn text-right'>21</td></tr><tr><td class='valn'>Max length</td><td class='valn text-right'>2 883</td></tr><tr><td class='valn'>Mean length</td><td class='valn text-right'>394.9</td></tr><tr><td class='valn'>Min depth</td><td class='valn text-right'>2</td></tr><tr><td class='valn'>Max depth</td><td class='valn text-right'>2 376</td></tr><tr><td class='valn'>Mean depth</td><td class='valn text-right'>46.91</td></tr></tbody></table>
            </div>
        </div>
    
        <div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-4">
            <h5>Top 10 of contigs with the highest number of detected introns</h5>
                <span class="anchor" id="too_complex_detected"></span>

            <table class="table table-striped table-bordered table-sm mb-0 " style="width:100%">
        	    <tbody>
<tr><td class='valn'>GJJD01002985.1</td><td class='valn text-right'>2 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01085655.1</td><td class='valn text-right'>2 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01000026.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01000401.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01000437.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01000474.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01000541.1</td><td class='valn text-right'>1 PASS; 2 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01002305.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01002766.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr><tr><td class='valn'>GJJD01003126.1</td><td class='valn text-right'>1 PASS; 0 DP; 0 IO; 0 SS; 0 LEN</td></tr></tbody></table>  
            </div>
        </div>

        <div class="d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center mt-5 pb-2 border-bottom">
            <h1 class="h4">Glossary</h1>
                <span class="anchor" id="glossary"></span>
        </div>
		<div class="d-flex">
            <div class="mt-4 mr-0 pl-0 col-md-12">
                <span class="anchor" id="glossary"></span>
                This glossary reflects terms used in this report:<br>
                <strong>Candidat</strong>        : Detected/predicted introns filtrered by canonical junction, length and depth.</br>
                <strong>Canonic junction</strong>: GT_AG or CT_AC for donor and acceptor sites.</br>
                <strong>Contig</strong>          : A set of overlapping DNA segments that together represent a consensus region of DNA (source :Wikipedia)</br>
                <strong>Exon</strong>            : Coding part of DNA that is translated into protein.</br>
                <strong>Features</strong>        : A sequenced chain of nucleic acids which provide protein(s). (source:https://biologydictionary.net/gene/)</br>
                <strong>Fragment</strong>        : During sequencing, the fragment is a link between fowarded and reverse reads.</br>                
                <strong>Intron</strong>          : Any nucleotide sequence removed by RNA splicing during maturation of RNA. (source : Wikipedia)</br>
                <strong>Reads</strong>           : An inferred sequence of base pairs corresponding to all or part of a single DNA fragment.  (source : Wikipedia)</br>
                <strong>Sequence</strong>        : Base pair probabilities of a DNA fragment.</br>
                <strong>Splice events</strong>   : a process where exons may be included within or excluded from the final mRNA.</br>
                <strong>Trimmed contig</strong>  : Contig without retained intron(s).</br>
                <br>
                <img src="https://forgemia.inra.fr/emilien.lasguignes/intronSeeker/-/raw/master/doc/IntronSeeker-glossary.png" alt="intronSeeker glossary" style="width:825;height:245;">
            </div>
        </div>  
<br/>
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