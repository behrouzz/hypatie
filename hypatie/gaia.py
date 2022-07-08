# https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access#CommandLine_DataLink

from urllib.request import urlopen, urlretrieve
import pandas as pd
import json
import requests



BASE = 'https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE='

def datalink_url(retrieval_type,
           source_id,
           #format='fits'
           data_structure=None):

    valid_rt = ['EPOCH_PHOTOMETRY', 'XP_SAMPLED', 'XP_CONTINUOUS',
                 'MCMC_GSPPHOT', 'MCMC_MSC',  'RVS', 'ALL']
    valid_ds = ['INDIVIDUAL','COMBINED','RAW']

    # Check source_id => working on it

    # Check retrieval_type
    if retrieval_type not in valid_rt:
        raise Exception(f'retrieval_type not valid! Options::\n{valid_rt}')

    # Check data_structue
    if data_structure is not None:
        if data_structure not in valid_ds:
            raise Exception(f'data_structure not valid! Options::\n{valid_ds}')
        else:
            data_structure = '&DATA_STRUCTURE=' + data_structure
    else:
        data_structure = ''

    url = BASE + retrieval_type+ '&' + f'ID=Gaia+DR3+{source_id}' + \
          data_structure + '&' + 'RELEASE=Gaia+DR3&FORMAT=fits'

    return url


def download(url, filename):
    urlretrieve(url, filename)
    

retrieval_type = 'EPOCH_PHOTOMETRY'
source_id = 30343944744320
data_structure = 'RAW'

url = datalink_url(retrieval_type, source_id)#, data_structure)
print(url)

"""
# https://dms.cosmos.esa.int/COSMOS/doc_fetch.php?id=2779219

# or: format(num, 'b')

def get_bin(num):
    num_int = int(bin(num)[2:])
    b = str(num_int).zfill(64)
    return b

def heal(num):
    b = get_bin(num)
    h = b[:29]
    b12 = int(b[:29],2)
    b11 = int(b[:27],2)
    b10 = int(b[:25],2)
    b09 = int(b[:23],2)
    b08 = int(b[:21],2)
    b07 = int(b[:19],2)
    b06 = int(b[:17],2)
    b05 = int(b[:15],2)
    b04 = int(b[:13],2)
    b03 = int(b[:11],2)
    b02 = int(b[:9],2)
    b01 = int(b[:7],2)
    b00 = int(b[:5],2)
    return [b00,b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12]

num = 196656661174768512

ls = heal(num)
b00,b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12 = ls

nsidS = [2**i for i in [0,1,2,3,4,5,6,7,8,9,10,11,12]]
npixS = [12*(nsid**2) for nsid in nsidS]

for ind, v in enumerate(ls):
    print(f'Pixel n° {v} from rang 0 to {npixS[ind]}')

# unique id
uS = []
for i in range(len(ls)):
    uS.append(ls[i] + 4*(nsidS[i]**2))

# https://stackoverflow.com/questions/29702010/healpy-pix2ang-convert-from-healpix-index-to-ra-dec-or-glong-glat


"""
