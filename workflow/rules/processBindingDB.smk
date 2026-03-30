import sys
from pathlib import Path
from damply import dirs
import pandas as pd
import requests
import zipfile


info = config['binding_db']
subset = info['subset']
version = info['version']

rule download_BindingDB:
	params:
		base_url = info['base_url']

	output:
		zipped_data = dirs.RAWDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_tsv.zip"

	run: 
		outpath = dirs.RAWDATA / "BINDING_DB"
		fname = f"BindingDB_{subset}_{version}_tsv.zip"
		
		url = "_".join([params.base_url,subset, version,"tsv.zip"])
		print(url)
		Path(outpath).mkdir(parents=True,exist_ok=True)
		print("about to fetch bdb")
		res = requests.get(url)
		
		with open(outpath / fname,"wb") as f:
			f.write(res.content)


rule process_BindingDB:
	params:
		organisms = info['processing']['keep_organisms'],
		org_col = info['processing']['organism_col'],
		useful_cols = info['processing']['keep_cols'],
		assay_col = info['processing']['assay_col'],
		cid_col = info['processing']['cid_col'],
		output_cols = info['processing']['output_cols']

	input: 
		raw_zip = dirs.RAWDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_tsv.zip",
			

	output:
		cleaned_data = dirs.PROCDATA / "BINDING_DB" / f"BindingDB_{subset}_{version}_cleaned.csv"

		
	run:
		Path(output.cleaned_data).parent.mkdir(parents=True, exist_ok=True)
		with zipfile.ZipFile(input.raw_zip,'r') as zf:
			zf.extractall(dirs.RAWDATA / "BINDING_DB")
		

		data = pd.read_csv(dirs.RAWDATA / "BINDING_DB" / "BindingDB_All.tsv",usecols = params.useful_cols,sep="\t")

		data = data[(data[params.org_col].isin(params.organisms)) & (data[params.assay_col].isnull())][params.output_cols]
		
		data.dropna(inplace=True, subset =[params.cid_col])
		
		data[params.cid_col] =[int(cid) for cid in data[params.cid_col].values]
		
		data.to_csv(output.cleaned_data,index=False)

