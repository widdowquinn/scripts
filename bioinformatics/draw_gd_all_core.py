#!/usr/bin/env python
#
# draw_gd_all_core.py
#
# (c) The James Hutton Institute 2013
# Author: Leighton Pritchard
#
# Contact:
# leighton.pritchard@hutton.ac.uk
#
# Leighton Pritchard,
# Information and Computing Sciences,
# James Hutton Institute,
# Errol Road,
# Invergowrie,
# Dundee,
# DD6 9LH,
# Scotland,
# UK
#
# The MIT License
#
# Copyright (c) 2010-2014 The James Hutton Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

""" Script to draw a simple pairwise GenomeDiagram figure with connectors,
    based on i-ADHoRe data
"""

# builtins
import os
from itertools import chain

# Biopython
from Bio.Graphics import GenomeDiagram as gd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Reportlab
from reportlab.lib.units import cm
from reportlab.lib import colors

# local
from ColorSpiral import get_color_dict, get_colors
from iadhore import IadhoreData

# Genome data: (id, GenBank file location)
# We list organisms in the order we want to present them, which is from
# 'outside-in' on the clade tree (GenomeDiagram orders from the bottom-up)
orgs = ('DPA2511', 'DPA703',
        'DXXRW240', 'DXXDW0440',
        'DZE3531', 'DZERW192', 'DZE586', 'DZE3532', 'DZEMK19', 'DZE2538',
        'DCH1591', 'DCH3533', 'DCH516', 'DCH402',
        'DXY569',
        'DSOIPO2222', 'DSOGBBC2040', 'DSOMK10', 'DSOMK16',
        'DXX3274', 'DXXMK7',
        'DDA898', 'DDA2976', 'DDA3937', 'DDA3537',
        'DDI3534', 'DDIIPO980', 'DDI453', 'DDIGBBC2039')

# List of organisms that need to be reverse complemented
reverse = ('DCH1591', 'DCH3533', 'DCH516', 'DCH402', 'DPA2511', 'DPA703')

gbkdir = '/mnt/synology_Dickeya_sequencing/Annotations/' +\
         '20120813_local_genbank_annotation/'
refdir = '/mnt/synology_Dickeya_sequencing/NCBI_GenBank_reference/'
genome_data = {'DPA2511': (gbkdir, 'NCPPB_2511_draft.gbk'),
               'DPA703': (refdir, 'NC_012880.gbk'),
               'DXXRW240': (gbkdir, 'CSL_RW240_draft.gbk'),
               'DXXDW0440': (gbkdir, 'DW_0440_draft.gbk'),
               'DZE3531': (gbkdir, 'NCPPB_3531_draft.gbk'),
               'DZERW192': (gbkdir, 'CSL_RW192_draft.gbk'),
               'DZE586': (refdir, 'NC_013592.gbk'),
               'DZE3532': (gbkdir, 'NCPPB_3532_draft.gbk'),
               'DZEMK19': (gbkdir, 'MK19_draft.gbk'),
               'DZE2538': (gbkdir, 'NCPPB_2538_draft.gbk'),
               'DCH1591': (refdir, 'NC_012912.gbk'),
               'DCH3533': (gbkdir, 'NCPPB_3533_draft.gbk'),
               'DCH516': (gbkdir, 'NCPPB_516_draft.gbk'),
               'DCH402': (gbkdir, 'NCPPB_402_draft.gbk'),
               'DXY569': (gbkdir, 'NCPPB_569_draft.gbk'),
               'DSOIPO2222': (gbkdir, 'IPO_2222_draft.gbk'),
               'DSOGBBC2040': (gbkdir, 'GBBC2040_draft.gbk'),
               'DSOMK10': (gbkdir, 'MK10_draft.gbk'),
               'DSOMK16': (gbkdir, 'MK16_draft.gbk'),
               'DXX3274': (gbkdir, 'NCPPB_3274_draft.gbk'),
               'DXXMK7': (gbkdir, 'MK7_draft.gbk'),
               'DDA898': (gbkdir, 'NCPPB_898_draft.gbk'),
               'DDA2976': (gbkdir, 'NCPPB_2976_draft.gbk'),
               'DDA3937': (refdir, 'NC_014500.gbk'),
               'DDA3537': (gbkdir, 'NCPPB_3537_draft.gbk'),
               'DDI3534': (gbkdir, 'NCPPB_3534_draft.gbk'),
               'DDIIPO980': (gbkdir, 'IPO_980_draft.gbk'),
               'DDIGBBC2039': (gbkdir, 'GBBC2039_draft.gbk'),
               'DDI453': (gbkdir, 'NCPPB_453_draft.gbk')
               }

# Create GenomeDiagram image
gdd = gd.Diagram("Dickeya core collinear regions", x=0.01, y=0.005)
tracks = {}
featuresets = {}
regionsets = {}
records = {}
track_level = 1
org_colours = get_color_dict(orgs, a=4)
for l, org in enumerate(orgs):
    # Load data
    filename = os.path.join(genome_data[org][0], genome_data[org][1])
    print "Loading %s" % filename
    records[org] = SeqIO.read(filename, 'genbank')
    if org in reverse:
        print "Reverse-complementing %s" % org
        records[org] = records[org].reverse_complement(annotations=True,
                                                       id=True, name=True,
                                                       description=True)
    # Set up tracks
    tracks[org] = gdd.new_track(2*l, name=org, greytrack=True,
                                greytrack_labels=10,
                                height=0.5,
                                start=0, end=len(records[org]))
    regionsets[org] = tracks[org].new_set(name="collinear regions")


# Convenience function for getting feature locations
def get_ft_loc(org, ft):
    for f in records[org].features:
        if f.type == 'CDS' and f.qualifiers['locus_tag'][0] == str(ft):
            return f.location.nofuzzy_start, f.location.nofuzzy_end

# Get data for crosslinks from i-ADHoRe results
data = IadhoreData(os.path.join('dickeya_all_output_params2',
                                'multiplicons.txt'),
                   os.path.join('dickeya_all_output_params2',
                                'segments.txt'))
full_leaves = data.get_multiplicons_at_level(29)
region_colours = list(get_colors(len(full_leaves), a=5, b=0.33,
                                 jitter=0.25))
for midx, m in enumerate(full_leaves):
    segments = data.get_multiplicon_segments(m)
    # Loop over the pairs of consecutive genomes in the table, and add
    # crosslinks for multiplicons
    print "Cross-linking multiplicon %d" % m
    for idx in range(1, len(orgs)):
        org1, org2 = orgs[idx-1], orgs[idx]
        org1loc = list(chain.from_iterable([get_ft_loc(org1, f) for f in
                                            segments[org1]]))
        org2loc = list(chain.from_iterable([get_ft_loc(org2, f) for f in
                                            segments[org2]]))
        org1ft = (tracks[org1], min(org1loc), max(org1loc))
        org2ft = (tracks[org2], min(org2loc), max(org2loc))
        # Need to create a colour rather than pass a tuple - unlike features.
        # Raise this bug in Biopython!
        c = colors.Color(region_colours[midx][0], region_colours[midx][1],
                         region_colours[midx][2])
        crosslink = gd.CrossLink(org1ft, org2ft, c)
        gdd.cross_track_links.append(crosslink)
        # Add feature to track (with transparency)
        # We add org1 here, then the final org when the looping's done
        f = SeqFeature(FeatureLocation(min(org1loc), max(org1loc)),
                       strand=None)
        regionsets[org1].add_feature(f, label=False, color=
                                     colors.Color(region_colours[midx][0],
                                                  region_colours[midx][1],
                                                  region_colours[midx][2],
                                                  0.5))
    # Finish off the cross-link features
    f = SeqFeature(FeatureLocation(min(org2loc), max(org2loc)),
                   strand=None)
    regionsets[org2].add_feature(f, label=False,
                                 color=colors.Color(region_colours[midx][0],
                                                    region_colours[midx][1],
                                                    region_colours[midx][2],
                                                    0.5))

# Add annotated organism features
for l, org in enumerate(orgs):
    print "Adding features for %s" % org
    featuresets[org] = tracks[org].new_set(name="CDS features")
    label_state = True
    for feature in [f for f in records[org].features if f.type == 'CDS']:
        label_state = not label_state
        featuresets[org].add_feature(feature, color=org_colours[org],
                                     label=False, sigil="ARROW",
                                     label_size=3)

# Render image
print "Rendering"
gdd.draw(format='linear', orientation='landscape',
         pagesize=(500*cm, (len(orgs)*91.4/29)*cm),
         fragments=1)
gdd.write('dickeya_core_collinear.pdf', 'PDF')
gdd.write('dickeya_core_collinear.png', 'PNG')
