import pandas as pd
import numpy as np
import argparse
import yaml
import logging

from logging.handlers import TimedRotatingFileHandler

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# create file handler which logs even debug messages
fh = TimedRotatingFileHandler('pathotype_identification.log', when='D', interval=1, backupCount=5)
fh.setLevel(logging.DEBUG)

# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)

# add the handlers to the logger
logger.addHandler(fh)

# File paths
#yaml_file = './designation.yml'
#input_file = './Monogenic_lines_10523.csv'
#output_file = './output_IRBL_Data.csv'

# Argument parser
parser = argparse.ArgumentParser(description='Identifying pathotype')
parser.add_argument('yaml_file', type=str, help='Path to YAML file')
parser.add_argument('input_file', type=str, help='Path to input file')
parser.add_argument('output_file', type=str, help='Path to output file')
args = parser.parse_args()

def convert_m_to_r_or_s(yaml_data, input_data):
    """
    Convert M to R or S
    """
    for key in yaml_data:
        if key in ['sh-B', 'ta2-Pi']:
            input_data.loc[input_data[key] == 'M', key] = 'R'
            logger.debug(f'Converted {key} M to R')
        else:
            input_data.loc[input_data[key] == 'M', key] = 'S'
            logger.debug(f'Converted {key} M to S')
    return input_data

def convert_sr_rs(yaml_data, input_data):
    """
    Convert S/R to S and R/S to R
    """
    for key in yaml_data:
        input_data.loc[input_data[key] == 'S/R', key] = 'S'
        logger.debug(f'Converted {key} S/R to S')
        input_data.loc[input_data[key] == 'R/S', key] = 'R'
        logger.debug(f'Converted {key} R/S to R')
        input_data.loc[input_data[key] == 'R/M', key] = 'R'
        logger.debug(f'Converted {key} R/M to R')
        input_data.loc[input_data[key] == 'R-M', key] = 'R'
        logger.debug(f'Converted {key} R-M to R')
        input_data.loc[input_data[key] == 'MS', key] = 'S'
        logger.debug(f'Converted {key} MS to S')
    return input_data

def process_data(yaml_file, input_file, output_file):
    try:
        # Load YAML file
        with open(args.yaml_file) as f:
            config = yaml.safe_load(f)

        # Extract keys from YAML file
        column_names = ['Designation'] + list(config.keys())

        # Load data from CSV file
        input_data = pd.read_csv(args.input_file)

        input_data = convert_m_to_r_or_s(yaml_file, input_data)
        input_data = convert_sr_rs(yaml_file, input_data)

        # Preprocess data
        prep_data = input_data.copy()
        prep_data['Designation'] = prep_data['Designation'].str.replace('IRBL', '')
        transposed = prep_data.T
        transposed.reset_index(inplace=True)
        transposed.columns = transposed.iloc[0]
        transposed.rename(columns=transposed.iloc[0], inplace=True)
        transposed.drop(index=transposed.index[0], axis=0, inplace=True)
        filtered = transposed[column_names]
        final_list = filtered.loc[:, ~filtered.columns.duplicated()].copy()

        # Convert M to R or S
        for key in config:
            if key in ['sh-B', 'ta2-Pi']:
                final_list.loc[final_list[key] == 'M', key] = 'R'
            else:
                final_list.loc[final_list[key] == 'M', key] = 'S'

        # Convert S/R to S and R/S to R
        for key in config:
            final_list.loc[final_list[key] == 'S/R', key] = 'S'
            final_list.loc[final_list[key] == 'R/S', key] = 'R'
            final_list.loc[final_list[key] == 'R/M', key] = 'R'
            final_list.loc[final_list[key] == 'R-M', key] = 'R'
            final_list.loc[final_list[key] == 'MS', key] = 'S'

        # Convert other cell values to NaN
        for key in config:
            final_list.loc[final_list[key] == '?', key] = np.nan
            final_list.loc[final_list[key] == '?R', key] = np.nan
            final_list.loc[final_list[key] == 'NaN', key] = np.nan
            final_list.loc[final_list[key] == '-', key] = np.nan
            final_list.loc[final_list[key] == 'nan', key] = np.nan
            final_list.loc[final_list[key] == 'D', key] = np.nan

        final_list.dropna(inplace=True)

        # Convert R and S to numeric value
        for key in config:
            final_list.loc[final_list[key] == 'S', key] = config[key]
            final_list.loc[final_list[key] == 'R', key] = 0

        # Add pathotype (eg. U73-i7-k177-z17-ta733)
        Ia = ['sh-S', 'b-B', 't-K59']
        #Ia = ['sh-B', 'b-B', 't-K59']
        Ib = ['LTH', 'a-A']
        II = ['i-F5', '3-CP4', '5-M']
        IIIa = ['ks-F5']
        IIIb = ['km-Ts', '1-CL', 'kh-K3']
        IIIc = ['k-Ka', 'kp-K60', '7-M']
        IVa = ['9-W']
        IVb = ['z-Fu', 'z5-CA', 'zt-T']
        Va = ['ta2-Pi', 'ta2-Re', '12-M']
        #Vb = ['ta-K1', 'ta-CP1']
        Vb = ['ta-K1', 'ta-CT2']
        Vc = ['19-A', '20-IR24']

        final_list['IA', 'IB', 'II', 'IIIA', 'IIIB', 'IIIC', 'IVA', 'IVB', 'VA', 'VB', 'VC'] = ''

        final_list['IA'] = final_list[Ia].sum(axis=1).astype(np.int64).astype(str)
        final_list['IB'] = final_list[Ib].sum(axis=1).astype(np.int64).astype(str)
        final_list['II'] = final_list[II].sum(axis=1).astype(np.int64).astype(str)
        final_list['IIIA'] = final_list[IIIa].sum(axis=1).astype(np.int64).astype(str)
        final_list['IIIB'] = final_list[IIIb].sum(axis=1).astype(np.int64).astype(str)
        final_list['IIIC'] = final_list[IIIc].sum(axis=1).astype(np.int64).astype(str)
        final_list['IVA'] = final_list[IVa].sum(axis=1).astype(np.int64).astype(str)
        final_list['IVB'] = final_list[IVb].sum(axis=1).astype(np.int64).astype(str)
        final_list['VA'] = final_list[Va].sum(axis=1).astype(np.int64).astype(str)
        final_list['VB'] = final_list[Vb].sum(axis=1).astype(np.int64).astype(str)
        final_list['VC'] = final_list[Vc].sum(axis=1).astype(np.int64).astype(str)

        final_list['Pathotype'] = 'U' + final_list['IA'] + final_list['IB'] + '-i' + \
                                final_list['II'] + '-k' + final_list['IIIA'] + \
                                final_list['IIIB'] + final_list['IIIC'] + '-z' + \
                                final_list['IVA'] + final_list['IVB'] + '-ta' + \
                                final_list['VA'] + final_list['VB'] + final_list['VC']

        grouped = final_list.sort_values(by=['Pathotype'])

        final_output = grouped[['Designation', 'Pathotype']].copy()

        # Save processed data to output file
        final_output.to_csv(args.output_file, index=False)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        logging.error(f"Error: {e}")
    except Exception as e:
        print(f"Error: {e}")
        logging.error(f"Error: {e}")