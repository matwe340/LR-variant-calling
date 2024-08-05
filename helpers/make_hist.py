from argparse import ArgumentParser
import os
import sys
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def main(input, out_folder, value='value', raw_data=None):
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
            
    if raw_data is None:
        for i, line in enumerate(input):
            line = line.strip()
            if i == 0:
                samples = [s for s in line.split('\t')]
                sample_hists = np.zeros((1000, len(samples)))
                continue
            
            for i, value in enumerate(line.split('\t')):
                if value == '.':
                    value = 0
                value = int(value)
                if value > sample_hists.shape[0]:
                    sample_hists = np.vstack((sample_hists, np.zeros((100, sample_hists.shape[1]))))
                
                sample_hists[value, i] += 1
        
        i = sample_hists.shape[0] - 1
        while i > 0 and np.sum(sample_hists[i, :]) == 0:
            i -= 1
        sample_hists = sample_hists[:i+1, :]
        max_value = i
        
        pd.DataFrame(sample_hists, columns=samples).to_csv(os.path.join(out_folder, 'raw_data.csv'), index=True)
    else:
        sample_hists = raw_data
        max_value = sample_hists.shape[0] - 1
        samples = sample_hists.columns
                        
    for i, sample in enumerate(samples):
        plt.bar(np.arange(max_value + 1), sample_hists[sample])
        plt.title(sample)
        plt.xlabel(value)
        plt.ylabel('Count')
        plt.savefig(os.path.join(out_folder, sample.split('/')[-1].split(':')[0] + '.png'))
        plt.close()


if __name__ == '__main__':
    parser = ArgumentParser()
    
    parser.add_argument('out_folder')
    parser.add_argument('--value', default='value')
    parser.add_argument('--input', default=None)
    parser.add_argument('--raw_data', default=None)
    
    args = parser.parse_args()
    
    input = None
    raw_data = None
    if args.input is None and args.raw_data is None:
        input = sys.stdin
    else:
        if args.input is not None:
            input = open(args.input)
        if args.raw_data is not None:
            raw_data = pd.read_csv(args.raw_data)
    
    main(input, args.out_folder, args.value, raw_data)
    
    if input is not None:
        input.close()