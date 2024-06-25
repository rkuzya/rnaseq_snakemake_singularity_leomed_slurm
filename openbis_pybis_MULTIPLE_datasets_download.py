### import packages ###
import pandas as pd
from pybis import Openbis
o = Openbis('https://openbis-sisrp-01.leomed.ethz.ch/openbis/webapp/eln-lims/?', verify_certificates=False)
#o = Openbis('example.com')          # https:// is assumed

import getpass
#password = getpass.getpass()
#o.login('admin', "changeit", save_token=True)
o.login('rkuzyakiv', "Lozynska2457", save_token=True)

### TO DOWNLOAD THE DATASETS of the certain type under the sample from the openBIS (into the Jupyter notebook) ###
##################################################################################################################
sample = o.get_sample('/OPENBIS_ONDEMAND/REGISTER_OPENBIS_ENTITIES_WITH_PYBIS/EXP1')

datasets = sample.get_datasets(type='RAW_DATA', start_with=0, count=10)
for dataset in datasets:
    print(dataset.props())
    print(dataset.file_list)
    dataset.download(
       destination = './input',        # download files to folder my_data/
	   create_default_folders = False, # ignore the /original/DEFAULT folders made by openBIS
	   wait_until_finished = False,    # download in background, continue immediately
	   workers = 10                    # 10 downloads parallel (default)
	)
dataset = datasets[0]



## GO TO openBIS to find the permID of the dataset you want to download

#ds = o.get_dataset("20231124145600817-63") # your PermID will be here
#ds = o.get_dataset("20240424150129475-101") # your PermID will be here
#ds.download(
#    destination = './input/', # download files to folder my_data/
#    create_default_folders = False, # ignore the /original/DEFAULT folders made by openBIS
#    wait_until_finished = False,    # download in background, continue immediately
#    workers = 10                    # 10 downloads parallel (default)
#)
#
#### the file is in the folder DOWNLOADED_DATA

