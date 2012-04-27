#!/usr/bin/env python -O

import xml.etree.ElementTree as eltree
import csv
import json

def measurement_dict(dataset_elem):
    meas_dict = dict()
    for elem in dataset_elem.findall('.//{http://phenome.jax.org}measurement'):
        id = int(elem.find('{http://phenome.jax.org}id').text[4:])
        title = elem.find('{http://phenome.jax.org}title').text
        abst = elem.find('{http://phenome.jax.org}abstract').text
        
        if id in meas_dict:
            raise Exception('found a duplicate measurement id: ' + id)
        
        meas_dict[id] = {'title': title, 'abstract': abst}
    
    return meas_dict

def dataset_dict(root_elem):
    dat_dict = dict()
    for elem in root_elem.findall('{http://phenome.jax.org}dataset'):
        id = int(elem.get('idnumber')[4:])
        title = elem.find('{http://phenome.jax.org}title').text
        meas = measurement_dict(elem)
        
        if id in dat_dict:
            raise Exception('found a duplicate dataset id: ' + id)
        
        if meas:
            dat_dict[id] = {'title': title, 'measurements': meas}
    
    return dat_dict

def read_pheno_dict(pheno_reader):
    """
    Read in the phenotype CSV file into a dictionary
    """
    
    # skips the header row
    pheno_reader.next()
    
    pheno_dict = dict()
    for row in pheno_reader:
        measnum = int(row[0])
        strain = row[1]
        sex = row[2]
        animal_id = row[3]
        value = None
        try:
            value = float(row[4])
        except:
            continue
        
        curr_pheno_row = {'animal_id': animal_id, 'strain': strain, 'sex': sex, 'value': value}
        
        if(measnum in pheno_dict):
            pheno_dict[measnum].append(curr_pheno_row)
        else:
            pheno_dict[measnum] = [curr_pheno_row]

    return pheno_dict

# main entry point for script
def main():
    datasets_root = eltree.parse('datasets_metadata.xml').getroot()
    pheno_reader = csv.reader(open('animaldatapoints.csv', 'rb'))
    print json.dumps(dataset_dict(datasets_root), indent=2)
    print json.dumps(read_pheno_dict(pheno_reader), indent=2)

if __name__ == '__main__':
    main()
