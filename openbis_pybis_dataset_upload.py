### import packages ###
import pandas as pd
from pybis import Openbis
o = Openbis('https://openbis-sisrp-01.leomed.ethz.ch/openbis/webapp/eln-lims/?', verify_certificates=False)
#o = Openbis('example.com')          # https:// is assumed

#import getpass
#password = getpass.getpass()
o.login('USERNAME', "PASSWORD", save_token=True)  

### REGISTER A NEW DATA SET (under the experimental steps) ###
##############################################################

# to get the experiment and sample use the information provided in the sample.
# see the value if the property Identifier (in the output of the previous cell) 
ds_new = o.new_dataset(
    type       = 'RAW_DATA',
    experiment = '/OPENBIS_ONDEMAND/REGISTER_OPENBIS_ENTITIES_WITH_PYBIS/DEMO_PYBIS', # ADAPT ACCORDING TO YOUR openBIS user training ID (number)
    sample     = '/OPENBIS_ONDEMAND/REGISTER_OPENBIS_ENTITIES_WITH_PYBIS/EXP1', # ADAPT ACCORDING TO YOUR openBIS user training ID (number) and EXPERIMENT
#    files      = ['./output/A_fastqc.html'],
    files      = ['./output'],
    props      = {'$name': 'fastqc report' }
)
ds_new.save()
