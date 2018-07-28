#!/usr/bin/env python
'''
  given tumour and normal vcf pairs, explore msi status
'''

import argparse
import logging
import sys

import cyvcf2
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
plt.rc('savefig', dpi=300)

import pandas as pd
import seaborn as sns

import intervaltree

rc={'font.size': 10, 'axes.labelsize': 10, 'legend.fontsize': 10.0, 'axes.titlesize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10}
sns.set(rc=rc)

#AFS = [0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
AFS = [0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5]
DPS = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400]
FIGSIZE = (20, 16)

def make_intervals(bed):
  intervals = {}
  total = 0
  for idx, line in enumerate(open(bed, 'r')):
    fields = line.strip('\n').split('\t')
    if len(fields) < 3:
      logging.warn('skipped line %i in %s', idx + 1, bed)
      continue
    chr, start, finish = fields[:3]
    if chr not in intervals:
      intervals[chr] = intervaltree.IntervalTree()
    if len(intervals[chr][int(start):int(finish)]) > 0:
      logging.info('overlap at %s %s %s', chr, start, finish)
    intervals[chr][int(start):int(finish)] = True
    total += int(finish) - int(start)
    if (idx + 1) % 100000 == 0:
      logging.debug('%i lines...', idx + 1)

  logging.info('filtering to %i bases', total)
  return (total, intervals)

def main(vcfs, samples, target, max_dp, snvs_only, bed, per_mb, min_agreement):
  logging.info('starting...')

  variants = {}
  ranges = {'max_dp': 0, 'min_dp': 1e9, 'max_af': 0, 'min_af': 1.0 }

  if bed is not None:
    size, intervals = make_intervals(bed)

  for sample, vcf in zip(samples, vcfs):
    logging.info('reading %s from %s...', sample, vcf)
    vcf_in = cyvcf2.VCF(vcf)
    variant_count = skipped_pass = skipped_snvs = skipped_bed = included = 0
    sample_id = vcf_in.samples.index(sample)
    for variant_count, variant in enumerate(vcf_in):
      # GL000220.1      135366  .       T       C       .       LowEVS;LowDepth SOMATIC;QSS=1;TQSS=1;NT=ref;QSS_NT=1;TQSS_NT=1;SGT=TT->TT;DP=2;MQ=60.00;MQ0=0;ReadPosRankSum=0.00;SNVSB=0.00;SomaticEVS=0.71    DP:FDP:SDP:SUBDP:AU:CU:GU:TU    1:0:0:0:0,0:0,0:0,0:1,1 1:0:0:0:0,0:1,1:0,0:0,0
      if (variant_count + 1 ) % 100000 == 0:
        logging.info('reading %s: %i variants processed...', vcf, variant_count + 1)

      if len(variant.ALT) > 1:
        logging.debug('%s: variant %i is multi-allelic', vcf, variant_count + 1)

      if bed is not None and (variant.CHROM not in intervals or len(intervals[variant.CHROM][variant.POS]) == 0):
        logging.debug('%s: skipped due to bed %s %s', vcf, variant.CHROM, variant.POS)
        skipped_bed += 1
        continue

      if snvs_only and len(variant.ALT[0]) != len(variant.REF): 
        logging.debug('%s: skipped indel %s -> %s', vcf, variant.REF, variant.ALT[0])
        skipped_snvs += 1
        continue

      if variant.FILTER is not None: # PASS only
        skipped_pass += 1
        continue

      dp = variant.INFO["DP"]
      af = variant.format("AF")[sample_id][0]

      ranges = {
        'max_dp': max(ranges['max_dp'], dp),
        'min_dp': min(ranges['min_dp'], dp),
        'max_af': max(ranges['max_af'], af),
        'min_af': min(ranges['min_af'], af)
      }

      if max_dp is not None and dp > max_dp:
        continue

      included += 1

      variant_id = '{}/{}'.format(variant.CHROM, variant.POS)
      if variant_id not in variants:
        variants[variant_id] = [dp, af, 1]
      else:
        existing_dp, existing_af, existing_count = variants[variant_id]
        variants[variant_id] = [min(dp, existing_dp), min(af, existing_af), existing_count + 1]

    logging.info('reading %s: processed %i variants. skipped pass %i. skipped indel %i. skipped bed %i. included %i', vcf, variant_count + 1, skipped_pass, skipped_snvs, skipped_bed, included)

  counts = {}
  totals = {}
  for af in AFS:
    for dp in DPS:
      counts['{}|{}'.format(af, dp)] = 0
      totals['{}|{}'.format(af, dp)] = 0

  for variant in variants:
    dp, af, count = variants[variant]
    for cell in counts:
      af_min, dp_min = [float(x) for x in cell.split('|')]
      if af > af_min and dp > dp_min:
        totals[cell] += 1
        if count >= min_agreement:
          counts[cell] += 1

  sys.stdout.write('\n'.join(['{}\t{}'.format(k, ranges[k]) for k in ranges]))

  #sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format('AF>', 'DP>', 'Intersect', 'Total', '%Intersect'))
  #for cell in sorted(counts):
  #  if totals[cell] > 0:
  #    concordance = counts[cell] / totals[cell]
  #  else:
  #    concordance = 0
  #  af_min, dp_min = cell.split('|')
  #  sys.stdout.write('{}\t{}\t{}\t{}\t{:.1f}\n'.format(af_min, dp_min, counts[cell], totals[cell], concordance * 100))

  # make df
  df = pd.DataFrame(columns=['DP'] + AFS)
  ann = pd.DataFrame(columns=AFS)
  for dp in DPS:
    heatmap = [0] * len(AFS)
    ann_row = [0] * len(AFS)
    for idx, af in enumerate(AFS):
      if totals['{}|{}'.format(af, dp)] == 0:
        heatmap[idx] = 0
        ann_row[idx] = '0'
      else:
        heatmap[idx] = counts['{}|{}'.format(af, dp)] / totals['{}|{}'.format(af, dp)] * 100 # determines colour
        if bed is not None and per_mb:
          ann_row[idx] = '{}%\n{:.1f}\n{:.1f}'.format(int(heatmap[idx]), 1000000 * counts['{}|{}'.format(af, dp)] / size, 1000000 * totals['{}|{}'.format(af, dp)] / size) # what's displayed
        else:
          ann_row[idx] = '{}%\n{}\n{}'.format(int(heatmap[idx]), counts['{}|{}'.format(af, dp)], totals['{}|{}'.format(af, dp)]) # what's displayed
    series = pd.Series(data=[dp] + heatmap, index=df.columns)
    ann_series = pd.Series(data=ann_row, index=ann.columns)
    df = df.append(series, ignore_index=True)
    ann = ann.append(ann_series, ignore_index=True)

  df.set_index('DP', inplace=True)
  logging.info('plotting...')
  plt.figure(figsize=FIGSIZE )
  plot = sns.heatmap(df, annot=ann, fmt='')
  if bed is not None and per_mb:
    plot.set(xlabel='Minimum AF', ylabel='Minimum DP', title='Concordance % and muts/Mb for {}'.format(', '.join(samples)))
  else:
    plot.set(xlabel='Minimum AF', ylabel='Minimum DP', title='Concordance % for {}'.format(', '.join(samples)))
  figure = plot.get_figure()
  figure.savefig(target)

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Assess MSI')
  parser.add_argument('--vcfs', required=True, nargs='+', help='vcfs')
  parser.add_argument('--samples', required=True, nargs='+', help='vcfs')
  parser.add_argument('--target', required=False, default='concordance.png', help='image')
  parser.add_argument('--max_dp', required=False, type=int, help='max dp')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  parser.add_argument('--snvs_only', action='store_true', help='exclude indels')
  parser.add_argument('--bed', required=False, help='limit to provided regions')
  parser.add_argument('--per_mb', action='store_true', help='show per mb (if bed file included)')
  parser.add_argument('--min_agreement', required=False, type=int, help='minimum number of vcfs that must agree (default is all)')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  if args.min_agreement is None:
    min_agreement = len(args.vcfs)
  else:
    min_agreement = args.min_agreement
  main(args.vcfs, args.samples, args.target, args.max_dp, args.snvs_only, args.bed, args.per_mb, min_agreement)
