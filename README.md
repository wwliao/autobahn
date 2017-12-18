# svpipe
SV detection pipeline. Can run Delly, Lumpyexpress, MELT, Sveltor, etc. 

Can be run on a single tumor/normal pair
```bash
python run_svtools.py -h
```

Can be used with LSF for batch submission
```bash
python submit_svJobs.py -h
```

svtools_config provides necessary path to various softwares used. 
