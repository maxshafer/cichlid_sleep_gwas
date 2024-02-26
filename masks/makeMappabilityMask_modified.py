#!/usr/bin/env python

import gzip
import sys

class MaskGenerator:
  def __init__(self, filename, chr):
    self.lastCalledPos = -1
    self.lastStartPos = -1
    sys.stderr.write("making mask {}\n".format(filename))
    self.file = gzip.open(filename, "w")
    self.chr = chr
  
  # assume 1-based coordinate, output in bed format
  def addCalledPosition(self, pos):
    if self.lastCalledPos == -1:
      self.lastCalledPos = pos
      self.lastStartPos = pos
    elif pos == self.lastCalledPos + 1:
      self.lastCalledPos = pos
    else:
      self.file.write("{}\t{}\t{}\n".format(self.chr, self.lastStartPos - 1, self.lastCalledPos))
      self.lastStartPos = pos
      self.lastCalledPos = pos

with open("/home/ayasha/projects/def-mshafer/genome/masks/all_chr_mask_l100_s1_r0.9.fa", "r") as f:
  for line in f:
    if line.startswith('>'):
      chr = line.split()[0][1:]
      mask = MaskGenerator("/home/ayasha/projects/def-mshafer/genome/masks/chr_{}_mask_l100_s1_r0.9.bed.gz".format(chr), chr)
      pos = 0
      continue
    for c in line.strip():
      pos += 1
      if pos % 1000000 == 0:
        sys.stderr.write("processing pos:{}\n".format(pos))
      if c != "3":
        mask.addCalledPosition(pos)
      
